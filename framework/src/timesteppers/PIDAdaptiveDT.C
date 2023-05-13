//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "PIDAdaptiveDT.h"
#include "Transient.h"

registerMooseObject("MooseApp", PIDAdaptiveDT);

InputParameters
PIDAdaptiveDT::validParams()
{
  InputParameters params = BaseAdaptiveDT::validParams();
  params.addClassDescription("Determine time steps based on the PID control theory.");
  params.addRequiredParam<Real>("tolerance", "Tolerance of relative error below which local error estimate is maintained by the PID control scheme");
  params.addRequiredParam<Real>("K_integral", "The coefficient multiplying the integral term");
  params.addRequiredParam<Real>("K_proportional", "The coefficient multiplying the difference term");
  params.addRequiredParam<Real>("K_derivative", "The coefficient multiplying the derivative term");
  params.addParam<bool>("check_aux_error", false, "Whether to check auxiliary system variables to estimate local error");
  params.addParam<std::vector<std::string>>("check_variables", "List of variables for which local error is checked");
  params.addParam<int>("local_error_order", 2, "Order of local error estimate");
  params.addParam<int>("start_adapting", 2, "when to start taking adaptive time steps");
  return params;
}

PIDAdaptiveDT::PIDAdaptiveDT(const InputParameters & parameters)
  : BaseAdaptiveDT(parameters),
    _tolerance(getParam<Real>("tolerance")),
    _Kint(getParam<Real>("K_integral")),
    _Kpro(getParam<Real>("K_proportional")),
    _Kder(getParam<Real>("K_derivative")),
    _start_adapting(getParam<int>("start_adapting")),
    _local_error_order(getParam<int>("local_error_order")),
    _nl(_fe_problem.getNonlinearSystemBase()),
    _aux(_fe_problem.getAuxiliarySystem()),
    _check_aux_error(getParam<bool>("check_aux_error")),
    _check_variables(getParam<std::vector<std::string>>("check_variables"))
{
  _local_error_estimate_current = 0.;
  _local_error_estimate_old = 0.;
  _local_error_estimate_older = 0.;
  _local_error_estimate_older_save = 0.;
  _R = 1.;
  _adaptive_timestepping = true;

}

Real
PIDAdaptiveDT::computeDT()
{
  Real dt = _dt_old;
  computeAdaptiveDT(dt);
  return dt;
}

void
PIDAdaptiveDT::computeAdaptiveDT(Real & dt, bool allowToGrow, bool allowToShrink)
{
  if (_R == std::numeric_limits<Real>::max())
  {
    dt *= _R;
    return;
  }

  dt *= std::pow(_R, _Kint);
  if (_local_error_estimate_old > 0. && _local_error_estimate_older > 0.)
  {
    Real r1 = _local_error_estimate_old / _local_error_estimate_current;
    Real r2 = r1 * _local_error_estimate_old / _local_error_estimate_older;
    dt *= std::pow(r1, _Kpro) * std::pow(r2, _Kder);
  }
}

void
PIDAdaptiveDT::step()
{
  TimeStepper::step();

  if (_converged)
  {
    if (_t_step >= _start_adapting)
    {
      _local_error_estimate_older_save = _local_error_estimate_older;
      _local_error_estimate_older = _local_error_estimate_old;
      _local_error_estimate_old = _local_error_estimate_current;
      _local_error_estimate_current = localErrorEstimate();
      if (_local_error_estimate_current == 0.)
        _R = std::numeric_limits<Real>::max();
      else
        _R = _tolerance / _local_error_estimate_current;
      _cutback_factor = std::max(0.5, _R);
      _console << "Local Error Estimate: " << _local_error_estimate_current << " Local Error Estimate/Tolerance " << 1. / _R << std::endl;
    }
    else
      _R = 1.;
  }
  else
    _R = 0.;
}

bool
PIDAdaptiveDT::converged() const
{
  if (!_converged)
    return false;
  if (_R >= 1.0)
    return true;
  else
    return false;
}

void
PIDAdaptiveDT::rejectStep()
{
  BaseAdaptiveDT::rejectStep();
  _local_error_estimate_current = _local_error_estimate_old;
  _local_error_estimate_old = _local_error_estimate_older;
  _local_error_estimate_older = _local_error_estimate_older_save;
}

void
PIDAdaptiveDT::acceptStep()
{
  BaseAdaptiveDT::acceptStep();
}

void
PIDAdaptiveDT::postSolve()
{
  BaseAdaptiveDT::postSolve();
  if (_R < 1.0)
    _console << "Marking last solve not converged " << _local_error_estimate_current << " > " << _tolerance << std::endl;
}

Real
PIDAdaptiveDT::localErrorEstimate()
{
  Real err = 0.;
  if (isParamValid("check_variables"))
  {
    for (auto varname : _check_variables)
    {
      auto & var = _fe_problem.getStandardVariable(0, varname);
      Real rt = 0.;
      if (_local_error_order == 2)
        rt = getCurrentDT() / _dt_old;
      Real sumdiff = 0.;
      Real sum = 0.;
      if (var.isNodal())
      {
        for (const auto & node : _fe_problem.mesh().getMesh().local_node_ptr_range())
        {
          if (node->n_dofs(var.sys().number(), var.number()) == 0)
            continue;

          Real x = var.getNodalValue(*node);
          Real diff = 0.;
          if (((_t_step <= 2) && (_local_error_order == 2)) || _local_error_order == 1)
            diff = x - var.getNodalValueOld(*node);
          else if (_local_error_order == 2)
            diff = x - var.getNodalValueOld(*node) * (1. + rt) + var.getNodalValueOlder(*node) * rt;
          else
            mooseError("Local error order greater than 2 has not been implemented yet.");
          sumdiff += diff * diff;
          sum += x * x;
        }
      }
      else
      {
        for (const auto & elem : _fe_problem.mesh().getMesh().active_local_element_ptr_range())
        {
          if (elem->n_dofs(var.sys().number(), var.number()) == 0)
            continue;

          Real x = var.getElementalValue(elem);
          Real diff = 0.;
          if (((_t_step <= 2) && (_local_error_order == 2)) || _local_error_order == 1)
            diff = x - var.getElementalValueOld(elem);
          else if (_local_error_order == 2)
            diff = (x - var.getElementalValueOld(elem) * (1. + rt) + var.getElementalValueOlder(elem) * rt) * 0.5;
          else
            mooseError("Local error order greater than 2 has not been implemented yet.");
          sumdiff += diff * diff;
          sum += x * x;
        }
      }
      _communicator.sum(sumdiff);
      _communicator.sum(sum);
      err = std::max(err, std::sqrt(sumdiff) / std::sqrt(sum));
    }
  }
  else
  {
    if (_check_aux_error)
    {
      NumericVector<Number> & solution = _aux.solution();
      if (((_t_step <= 2) && (_local_error_order == 2)) || _local_error_order == 1)
      {
        NumericVector<Number> & sd = solution;
        sd -= _aux.solutionOld();
        err = sd.l2_norm() / solution.l2_norm();
      }
      else if (_local_error_order == 2)
      {
        Real rt = getCurrentDT() / _dt_old;
        NumericVector<Number> & temp = _aux.solutionOlder();
        temp *= rt;
        rt += 1.;
        rt *= -1.;
        NumericVector<Number> & sd = _aux.solutionOld();
        sd *= rt;
        sd += solution;
        sd += temp;
        sd *= 0.5;
        err = sd.l2_norm() / solution.l2_norm();
      }
      else
        mooseError("Local error order greater than 2 has not been implemented yet.");
    }
    else
    {
      const NumericVector<Number> & solution = *_nl.currentSolution();
      if (((_t_step <= 2) && (_local_error_order == 2)) || _local_error_order == 1)
      {
        NumericVector<Number> & sd = _nl.addVector("sd", false, PARALLEL);
        sd = solution;
        sd -= _nl.solutionOld();
        err = sd.l2_norm() / solution.l2_norm();
      }
      else if (_local_error_order == 2)
      {
        Real rt = getCurrentDT() / _dt_old;
        NumericVector<Number> & temp = _nl.solutionOlder();
        temp *= rt;
        rt += 1.;
        rt *= -1.;
        NumericVector<Number> & sd = _nl.solutionOld();
        sd *= rt;
        sd += solution;
        sd += temp;
        sd *= 0.5;
        err = sd.l2_norm() / solution.l2_norm();
      }
      else
        mooseError("Local error order greater than 2 has not been implemented yet.");
    }
  }
  return err;
}

// Real
// PIDAdaptiveDT::localErrorEstimate1(const NumericVector<Number> & s1, NumericVector<Number> & s2)
// {
//   NumericVector<Number> & sd = s1;
//   sd = s1;
//   sd -= s2;
//   return sd.l2_norm() / s1.l2_norm();
// }

// Real
// PIDAdaptiveDT::localErrorEstimate2(const NumericVector<Number> & s1, NumericVector<Number> & s2, NumericVector<Number> & s3)
// {
//   Real rt = getCurrentDT() / _dt_old;
//   NumericVector<Number> & temp = s3;
//   temp *= rt;
//   rt += 1.;
//   rt *= -1.;
//   NumericVector<Number> & sd = s2;
//   sd *= rt;
//   sd += s1;
//   sd += temp;
//   sd *= 0.5;
//   return sd.l2_norm() / s1.l2_norm();
// }
