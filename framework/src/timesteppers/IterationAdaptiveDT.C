//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "IterationAdaptiveDT.h"
#include "Transient.h"

registerMooseObject("MooseApp", IterationAdaptiveDT);

InputParameters
IterationAdaptiveDT::validParams()
{
  InputParameters params = BaseAdaptiveDT::validParams();
  params.addClassDescription("Adjust the timestep based on the number of iterations");
  params.addParam<int>("optimal_iterations",
                       "The target number of nonlinear iterations for adaptive timestepping");
  params.addParam<int>("iteration_window",
                       "Attempt to grow/shrink timestep if the iteration count "
                       "is below/above 'optimal_iterations plus/minus "
                       "iteration_window' (default = optimal_iterations/5).");
  params.addParam<unsigned>("linear_iteration_ratio",
                            "The ratio of linear to nonlinear iterations "
                            "to determine target linear iterations and "
                            "window for adaptive timestepping (default = "
                            "25)");
  params.addParam<bool>("reject_large_step",
                        false,
                        "If 'true', time steps that are too large compared to the "
                        "ideal time step will be rejected and repeated");
  params.addRangeCheckedParam<Real>("reject_large_step_threshold",
                                    0.1,
                                    "reject_large_step_threshold > 0 "
                                    "& reject_large_step_threshold < 1",
                                    "Ratio between the the ideal time step size and the "
                                    "current time step size below which a time step will "
                                    "be rejected if 'reject_large_step' is 'true'");
  return params;
}

IterationAdaptiveDT::IterationAdaptiveDT(const InputParameters & parameters)
  : BaseAdaptiveDT(parameters),
    _linear_iteration_ratio(isParamValid("linear_iteration_ratio")
                                ? getParam<unsigned>("linear_iteration_ratio")
                                : 25), // Default to 25
    _nl_its(declareRestartableData<unsigned int>("nl_its", 0)),
    _l_its(declareRestartableData<unsigned int>("l_its", 0)),
    _reject_large_step(getParam<bool>("reject_large_step")),
    _large_step_rejection_threshold(getParam<Real>("reject_large_step_threshold"))
{
  if (isParamValid("optimal_iterations"))
  {
    _adaptive_timestepping = true;
    _optimal_iterations = getParam<int>("optimal_iterations");

    if (isParamValid("iteration_window"))
      _iteration_window = getParam<int>("iteration_window");
    else
      _iteration_window = ceil(_optimal_iterations / 5.0);
  }
  else
  {
    if (isParamValid("iteration_window"))
      mooseError("'optimal_iterations' must be used for 'iteration_window' to be used");
    if (isParamValid("linear_iteration_ratio"))
      mooseError("'optimal_iterations' must be used for 'linear_iteration_ratio' to be used");
  }
}

Real
IterationAdaptiveDT::computeDT()
{
  Real dt = _dt_old;

  if (_cutback_occurred)
  {
    _cutback_occurred = false;
    if (_adaptive_timestepping)
    {
      // Don't allow it to grow this step, but shrink if needed
      bool allowToGrow = false;
      computeAdaptiveDT(dt, allowToGrow);
    }
  }
  else
    dt = BaseAdaptiveDT::computeDT();

  return dt;
}

bool
IterationAdaptiveDT::converged() const
{
  if (!_reject_large_step)
    return TimeStepper::converged();

  // the solver has not converged
  if (!TimeStepper::converged())
    return false;

  // we are already at dt_min or at the start of the simulation
  // in which case we can move on to the next step
  if (_dt == _dt_min || _t_step < 2)
    return true;

  // we get what the next time step should be
  Real dt_test = _dt;
  limitDTToPostprocessorValue(dt_test);

  // we cannot constrain the time step any further
  if (dt_test == 0)
    return true;

  // if the time step is much smaller than the current time step
  // we need to repeat the current iteration with a smaller time step
  if (dt_test < _dt * _large_step_rejection_threshold)
    return false;

  // otherwise we move one
  return true;
}

void
IterationAdaptiveDT::computeAdaptiveDT(Real & dt, bool allowToGrow, bool allowToShrink)
{
  const unsigned int growth_nl_its(
      _optimal_iterations > _iteration_window ? _optimal_iterations - _iteration_window : 0);
  const unsigned int shrink_nl_its(_optimal_iterations + _iteration_window);
  const unsigned int growth_l_its(_optimal_iterations > _iteration_window
                                      ? _linear_iteration_ratio *
                                            (_optimal_iterations - _iteration_window)
                                      : 0);
  const unsigned int shrink_l_its(_linear_iteration_ratio *
                                  (_optimal_iterations + _iteration_window));

  if (allowToGrow && (_nl_its < growth_nl_its && _l_its < growth_l_its))
  {
    // Grow the timestep
    dt *= _growth_factor;

    if (_verbose)
      _console << "Growing dt: nl its = " << _nl_its << " < " << growth_nl_its
               << " && lin its = " << _l_its << " < " << growth_l_its << " old dt: " << std::setw(9)
               << _dt_old << " new dt: " << std::setw(9) << dt << '\n';
  }
  else if (allowToShrink && (_nl_its > shrink_nl_its || _l_its > shrink_l_its))
  {
    // Shrink the timestep
    dt *= _cutback_factor;

    if (_verbose)
      _console << "Shrinking dt: nl its = " << _nl_its << " > " << shrink_nl_its
               << " || lin its = " << _l_its << " > " << shrink_l_its << " old dt: " << std::setw(9)
               << _dt_old << " new dt: " << std::setw(9) << dt << '\n';
  }

  _console << std::flush;
}

void
IterationAdaptiveDT::rejectStep()
{
  BaseAdaptiveDT::rejectStep();
}

void
IterationAdaptiveDT::acceptStep()
{
  BaseAdaptiveDT::acceptStep();

  _nl_its = _fe_problem.getNonlinearSystemBase().nNonlinearIterations();
  _l_its = _fe_problem.getNonlinearSystemBase().nLinearIterations();

  if ((_at_function_point || _executioner.atSyncPoint()) &&
      _dt + _timestep_tolerance < _executioner.unconstrainedDT())
  {
    _dt_old = _fe_problem.dtOld();
    _sync_last_step = true;

    if (_verbose)
      _console << "Sync point hit in current step, using previous dt for old dt: " << std::setw(9)
               << _dt_old << std::endl;
  }
}
