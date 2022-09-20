//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeStepper.h"
#include "LinearInterpolation.h"
#include "PostprocessorInterface.h"

class Function;
class PiecewiseBase;
class PiecewiseLinear;

/**
 * Adjust the timestep based on some control algorithm.
 * This time stepper allows the user to specify a function that limits the maximal time step change.
 * This time stepper allows the user to specify a limiting time step length through a postprocessor.
 */
class BaseAdaptiveDT : public TimeStepper, public PostprocessorInterface
{
public:
  static InputParameters validParams();

  BaseAdaptiveDT(const InputParameters & parameters);

  virtual void init() override;
  virtual void preExecute() override;

  virtual void rejectStep() override;
  virtual void acceptStep() override;

  virtual bool constrainStep(Real & dt) override;
  virtual void postSolve() override;

protected:
  virtual Real computeInitialDT() override;
  virtual Real computeDT() override;
  virtual Real computeFailedDT() override;

  virtual void computeAdaptiveDT(Real & dt) = 0;
  Real computeInterpolationDT();
  void limitDTByFunction(Real & limitedDT);
  void limitDTToPostprocessorValue(Real & limitedDT) const;

  Real & _dt_old;

  /// The dt from the input file.
  const Real _input_dt;

  bool & _tfunc_last_step;
  bool & _sync_last_step;

  /// adaptive timestepping is active if the optimal_iterations input parameter is specified
  bool _adaptive_timestepping;

  /// if specified, the postprocessor values used to determine an upper limit for the time step length
  std::vector<const PostprocessorValue *> _pps_value;

  std::vector<const Function *> _timestep_limiting_functions;
  std::vector<const PiecewiseBase *> _piecewise_timestep_limiting_functions;
  std::vector<const PiecewiseLinear *> _piecewise_linear_timestep_limiting_functions;

  /// time point defined in the piecewise function
  std::vector<Real> _times;

  Real _max_function_change;
  /// insert sync points at the time nodes of the _piecewise_timestep_limiting_function
  const bool _force_step_every_function_point;
  /// Set timestep size if previous timestep is synced with function
  const Real _post_function_sync_dt;

  std::set<Real> _tfunc_times;

  /// PiecewiseBase linear definition of time stepping
  LinearInterpolation _time_ipol;
  /// true if we want to use piecewise-defined time stepping
  const bool _use_time_ipol;

  /// grow the timestep by this factor
  const Real & _growth_factor;
  /// cut the timestep by this factor
  const Real & _cutback_factor;

  bool & _cutback_occurred;
  bool _at_function_point;

};
