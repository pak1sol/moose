//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "BaseAdaptiveDT.h"
#include "Function.h"
#include "PiecewiseLinear.h"
#include "Transient.h"
#include "NonlinearSystem.h"

#include <limits>
#include <set>

InputParameters
BaseAdaptiveDT::validParams()
{
  InputParameters params = TimeStepper::validParams();
  params.addClassDescription("Adjust the timestep based on time step control algorithm");
  params.addParam<std::vector<PostprocessorName>>("timestep_limiting_postprocessor",
                                                  "If specified, a list of postprocessor values "
                                                  "used as an upper limit for the "
                                                  "current time step length");
  params.addParam<std::vector<FunctionName>>(
      "timestep_limiting_function",
      "A list of 'PiecewiseBase' type functions used to control the timestep by "
      "limiting the change in the function over a timestep");
  params.addParam<Real>(
      "max_function_change",
      "The absolute value of the maximum change in timestep_limiting_function over a timestep");
  params.addParam<bool>("force_step_every_function_point",
                        false,
                        "Forces the timestepper to take "
                        "a step that is consistent with "
                        "points defined in the function");
  params.addRangeCheckedParam<Real>(
      "post_function_sync_dt",
      "post_function_sync_dt>0",
      "Timestep to apply after time sync with function point. To be used in "
      "conjunction with 'force_step_every_function_point'.");
  params.addRequiredParam<Real>("dt", "The default timestep size between solves");
  params.addParam<std::vector<Real>>("time_t", "The values of t");
  params.addParam<std::vector<Real>>("time_dt", "The values of dt");
  params.addParam<Real>("growth_factor",
                        2.0,
                        "Factor to apply to timestep if easy convergence (if "
                        "'optimal_iterations' is specified) or if recovering "
                        "from failed solve");
  params.addParam<Real>("cutback_factor",
                        0.5,
                        "Factor to apply to timestep if difficult "
                        "convergence (if 'optimal_iterations' is specified) "
                        "or if solution failed");
  params.declareControllable("growth_factor cutback_factor");

  return params;
}

BaseAdaptiveDT::BaseAdaptiveDT(const InputParameters & parameters)
  : TimeStepper(parameters),
    PostprocessorInterface(this),
    _dt_old(declareRestartableData<Real>("dt_old", 0.0)),
    _input_dt(getParam<Real>("dt")),
    _tfunc_last_step(declareRestartableData<bool>("tfunc_last_step", false)),
    _sync_last_step(declareRestartableData<bool>("sync_last_step", false)),
    _adaptive_timestepping(false),
    _pps_value(
        parameters.get<std::vector<PostprocessorName>>("timestep_limiting_postprocessor").size()),
    _timestep_limiting_functions(),
    _piecewise_timestep_limiting_functions(),
    _piecewise_linear_timestep_limiting_functions(),
    _times(0),
    _max_function_change(-1),
    _force_step_every_function_point(getParam<bool>("force_step_every_function_point")),
    _post_function_sync_dt(isParamValid("force_step_every_function_point") &&
                                   isParamValid("post_function_sync_dt")
                               ? getParam<Real>("post_function_sync_dt")
                               : 0.0),
    _tfunc_times(getParam<std::vector<Real>>("time_t").begin(),
                 getParam<std::vector<Real>>("time_t").end()),
    _time_ipol(getParam<std::vector<Real>>("time_t"), getParam<std::vector<Real>>("time_dt")),
    _use_time_ipol(_time_ipol.getSampleSize() > 0),
    _growth_factor(getParam<Real>("growth_factor")),
    _cutback_factor(getParam<Real>("cutback_factor")),
    _cutback_occurred(declareRestartableData<bool>("cutback_occurred", false)),
    _at_function_point(false)
{
  auto timestep_limiting_postprocessor_names =
      parameters.get<std::vector<PostprocessorName>>("timestep_limiting_postprocessor");
  for (size_t i = 0; i < _pps_value.size(); ++i)
    _pps_value[i] = &getPostprocessorValueByName(timestep_limiting_postprocessor_names[i]);

  if (isParamValid("timestep_limiting_function"))
    _max_function_change =
        isParamValid("max_function_change") ? getParam<Real>("max_function_change") : -1;
  else
  {
    if (isParamValid("max_function_change"))
      mooseError("'timestep_limiting_function' must be used for 'max_function_change' to be used");
    if (_force_step_every_function_point)
      mooseError("'timestep_limiting_function' must be used for 'force_step_every_function_point' "
                 "to be used");
  }

  if (!isParamValid("force_step_every_function_point") && isParamValid("post_function_sync_dt"))
    paramError("post_function_sync_dt",
               "Not applicable if 'force_step_every_function_point = false'");
}

void
BaseAdaptiveDT::init()
{
  if (isParamValid("timestep_limiting_function"))
  {
    std::set<Real> times;

    const auto tid = isParamValid("_tid") ? getParam<THREAD_ID>("_tid") : 0;
    for (const auto & name : getParam<std::vector<FunctionName>>("timestep_limiting_function"))
    {
      const auto * func = &_fe_problem.getFunction(name, tid);
      _timestep_limiting_functions.push_back(func);

      const auto * pfunc = dynamic_cast<const PiecewiseBase *>(func);
      _piecewise_timestep_limiting_functions.push_back(pfunc);

      if (pfunc)
      {
        const auto * plfunc = dynamic_cast<const PiecewiseLinear *>(pfunc);
        _piecewise_linear_timestep_limiting_functions.push_back(plfunc);

        const auto ntimes = pfunc->functionSize();
        for (unsigned int i = 0; i < ntimes; ++i)
          times.insert(pfunc->domain(i));
      }
      else
        mooseError("timestep_limiting_function must be a PiecewiseBase function");
    }
    _times.resize(times.size());
    std::copy(times.begin(), times.end(), _times.begin());

    mooseAssert(_timestep_limiting_functions.size() ==
                    _piecewise_timestep_limiting_functions.size(),
                "Timestep limiting function count inconsistency");
    mooseAssert(_piecewise_timestep_limiting_functions.size() ==
                    _piecewise_linear_timestep_limiting_functions.size(),
                "Timestep limiting function count inconsistency");
  }
}

void
BaseAdaptiveDT::preExecute()
{
  TimeStepper::preExecute();

  // Delete all tfunc times that are at or before the begin time
  while (!_tfunc_times.empty() && _time + _timestep_tolerance >= *_tfunc_times.begin())
    _tfunc_times.erase(_tfunc_times.begin());
}

Real
BaseAdaptiveDT::computeInitialDT()
{
  Real dt;
  if (_tfunc_last_step)
  {
    dt = _time_ipol.sample(_time_old);
    if (_verbose)
      _console << "Setting initial dt to value specified by dt function: " << std::setw(9) << dt
               << std::endl;
  }
  else
  {
    dt = _input_dt;
    if (_verbose)
      _console << "Setting initial dt to input value: " << std::setw(9) << dt << std::endl;
  }
  return dt;
}

Real
BaseAdaptiveDT::computeDT()
{
  Real dt = _dt_old;

  if (_tfunc_last_step)
  {
    _sync_last_step = false;
    dt = _time_ipol.sample(_time_old);

    if (_verbose)
      _console << "Setting dt to value specified by dt function: " << std::setw(9) << dt
               << std::endl;
  }
  else if (_sync_last_step)
  {
    _sync_last_step = false;
    if (_post_function_sync_dt)
      dt = _post_function_sync_dt;
    else
      dt = _dt_old;

    if (_verbose)
      _console << "Setting dt to value used before sync: " << std::setw(9) << dt << std::endl;
  }
  else if (_adaptive_timestepping)
    computeAdaptiveDT(dt);
  else if (_use_time_ipol)
    dt = computeInterpolationDT();
  else
  {
    dt *= _growth_factor;
    if (dt > _dt_old * _growth_factor)
      dt = _dt_old * _growth_factor;
    if (_verbose)
      _console << "Growing dt based on growth factor (" << _growth_factor
               << ") and previous dt before sync (" << _dt_old << ") : " << std::setw(9) << dt
               << std::endl;
  }

  return dt;
}

bool
BaseAdaptiveDT::constrainStep(Real & dt)
{
  bool at_sync_point = TimeStepper::constrainStep(dt);

  // Limit the timestep to postprocessor value
  limitDTToPostprocessorValue(dt);

  // Limit the timestep to limit change in the function
  limitDTByFunction(dt);

  // Adjust to the next tfunc time if needed
  if (!_tfunc_times.empty() && _time + dt + _timestep_tolerance >= *_tfunc_times.begin())
  {
    dt = *_tfunc_times.begin() - _time;

    if (_verbose)
      _console << "Limiting dt to sync with dt function time: " << std::setw(9)
               << *_tfunc_times.begin() << " dt: " << std::setw(9) << dt << std::endl;
  }

  return at_sync_point;
}

void
BaseAdaptiveDT::postSolve()
{
  bool tfunc_last_step = _tfunc_last_step;
  bool sync_last_step = _sync_last_step;
  bool cutback_occurred = _cutback_occurred;

  TimeStepper::postSolve();

  _tfunc_last_step = tfunc_last_step;
  _sync_last_step = sync_last_step;
  _cutback_occurred = cutback_occurred;
}

Real
BaseAdaptiveDT::computeFailedDT()
{
  _cutback_occurred = true;

  // Can't cut back any more
  if (_dt <= _dt_min)
    mooseError("Solve failed and timestep already at dtmin, cannot continue!");

  if (_verbose)
  {
    _console << "\nSolve failed with dt: " << std::setw(9) << _dt
             << "\nRetrying with reduced dt: " << std::setw(9) << _dt * _cutback_factor_at_failure
             << std::endl;
  }
  else
    _console << "\nSolve failed, cutting timestep." << std::endl;

  return _dt * _cutback_factor_at_failure;
}

void
BaseAdaptiveDT::limitDTToPostprocessorValue(Real & limitedDT) const
{
  if (_pps_value.size() != 0 && _t_step > 1)
  {
    Real limiting_pps_value = *_pps_value[0];
    unsigned int i_min = 0;
    for (size_t i = 1; i < _pps_value.size(); ++i)
      if (*_pps_value[i] < limiting_pps_value)
      {
        limiting_pps_value = *_pps_value[i];
        i_min = i;
      }

    if (limitedDT > limiting_pps_value)
    {
      if (limiting_pps_value < 0)
        mooseWarning(
            "Negative timestep limiting postprocessor '" +
            getParam<std::vector<PostprocessorName>>("timestep_limiting_postprocessor")[i_min] +
            "': " + std::to_string(limiting_pps_value));
      limitedDT = std::max(_dt_min, limiting_pps_value);

      if (_verbose)
        _console << "Limiting dt to postprocessor value. dt = " << limitedDT << std::endl;
    }
  }
}

void
BaseAdaptiveDT::limitDTByFunction(Real & limitedDT)
{
  Real orig_dt = limitedDT;
  const auto nfunc = _timestep_limiting_functions.size();
  Real restricted_step = std::numeric_limits<Real>::max();

  for (unsigned int j = 0; j < nfunc; ++j)
  {
    // Limit by function change for piecewise linear functions.
    if (_piecewise_linear_timestep_limiting_functions[j] && _max_function_change > 0)
    {
      const auto current_function_value =
          _piecewise_linear_timestep_limiting_functions[j]->value(_time, {});

      const auto ntimes = _piecewise_linear_timestep_limiting_functions[j]->functionSize();
      for (std::size_t next_time_index = 1; next_time_index < ntimes; ++next_time_index)
      {
        const auto next_time =
            _piecewise_linear_timestep_limiting_functions[j]->domain(next_time_index);

        // Skip ahead to find time point that is just past the current time.
        if (next_time + _timestep_tolerance <= _time)
          continue;

        // Find out how far we can go without exceeding the max function change.
        const auto next_function_value =
            _piecewise_linear_timestep_limiting_functions[j]->range(next_time_index);
        const auto change = std::abs(next_function_value - current_function_value);
        if (change > _max_function_change)
        {
          // Interpolate to find step.
          restricted_step =
              std::min(restricted_step, (_max_function_change / change) * (next_time - _time));
          break;
        }

        // Don't keep going if we've already passed the current limited step.
        if (next_time > _time + limitedDT)
          break;
      }
    }
    else if (_timestep_limiting_functions[j])
    {
      const Real old_value = _timestep_limiting_functions[j]->value(_time_old);
      Real new_value = _timestep_limiting_functions[j]->value(_time_old + limitedDT);
      Real change = std::abs(new_value - old_value);

      if (_max_function_change > 0.0 && change > _max_function_change)
        do
        {
          limitedDT /= 2.0;
          new_value = _timestep_limiting_functions[j]->value(_time_old + limitedDT);
          change = std::abs(new_value - old_value);
        } while (change > _max_function_change);
    }
  }

  if (restricted_step < limitedDT)
    limitedDT = std::max(_dt_min, restricted_step);

  _at_function_point = false;
  if (_force_step_every_function_point)
    for (unsigned int i = 0; i + 1 < _times.size(); ++i)
      if (_time >= _times[i] && _time < _times[i + 1])
      {
        if (limitedDT > _times[i + 1] - _time - _timestep_tolerance)
        {
          limitedDT = _times[i + 1] - _time;
          _at_function_point = true;
        }
        break;
      }

  if (_verbose && limitedDT != orig_dt)
  {
    if (_at_function_point)
      _console << "Limiting dt to match function point. dt = ";
    else
      _console << "Limiting dt to limit change in function. dt = ";

    _console << limitedDT << std::endl;
  }
}

Real
BaseAdaptiveDT::computeInterpolationDT()
{
  Real dt = _time_ipol.sample(_time_old);

  if (dt > _dt_old * _growth_factor)
  {
    dt = _dt_old * _growth_factor;

    if (_verbose)
      _console << "Growing dt to recover from cutback. "
               << " old dt: " << std::setw(9) << _dt_old << " new dt: " << std::setw(9) << dt
               << std::endl;
  }

  return dt;
}

void
BaseAdaptiveDT::rejectStep()
{
  TimeStepper::rejectStep();
}

void
BaseAdaptiveDT::acceptStep()
{
  TimeStepper::acceptStep();

  _tfunc_last_step = false;
  while (!_tfunc_times.empty() && _time + _timestep_tolerance >= *_tfunc_times.begin())
  {
    if (std::abs(_time - *_tfunc_times.begin()) <= _timestep_tolerance)
      _tfunc_last_step = true;

    _tfunc_times.erase(_tfunc_times.begin());
  }

  _dt_old = _dt;
}
