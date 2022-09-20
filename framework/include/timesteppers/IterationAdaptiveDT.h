//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "BaseAdaptiveDT.h"

/**
 * Adjust the timestep based on the number of iterations.
 * The user can specitfy an optimal_iterations number of non-linear iterations and an
 * iteration_window.
 * The time stepper attempts to increase the time step if the non-linear iteration count is below
 * optimal_iterations - iteration_window and it attempts to reduce the time step if the non-linear
 * iteration count is above optimal_iterations + iteration_window. Similar rules apply to the number
 * of linear iterations with the multiplier linear_iteration_ratio applied to optimal_iterations and
 * iteration_window.
 */
class IterationAdaptiveDT : public BaseAdaptiveDT
{
public:
  static InputParameters validParams();

  IterationAdaptiveDT(const InputParameters & parameters);

  virtual void rejectStep() override;
  virtual void acceptStep() override;

protected:
  virtual Real computeDT() override;
  virtual bool converged() const override;

  virtual void computeAdaptiveDT(Real & dt) override;

  /// Adapt the timestep to maintain this non-linear iteration count...
  int _optimal_iterations;
  /// ...plus/minus this value.
  int _iteration_window;
  /// use _optimal_iterations and _iteration_window multiplied with this factor for linear iterations
  const int _linear_iteration_ratio;

  /// Number of nonlinear iterations in previous solve
  unsigned int & _nl_its;
  /// Number of linear iterations in previous solve
  unsigned int & _l_its;

  /// Indicates whether we need to reject a time step much larger than its ideal size
  bool _reject_large_step;
  /// Threshold used to detect whether we need to reject a step
  double _large_step_rejection_threshold;

private:
  bool _allowToGrow;
  bool _allowToShrink = true;

};
