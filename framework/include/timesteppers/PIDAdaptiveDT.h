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
#include "NonlinearSystemBase.h"
#include "AuxiliarySystem.h"

/**
 * Adjust the timestep based on PID control theory.
 * The user can specitfy 
 */
class PIDAdaptiveDT : public BaseAdaptiveDT
{
public:
  static InputParameters validParams();

  PIDAdaptiveDT(const InputParameters & parameters);

  virtual void step() override;
  virtual bool converged() const override;
  virtual void rejectStep() override;
  virtual void acceptStep() override;
  virtual void postSolve() override;

protected:
  virtual Real computeDT() override;

  virtual void computeAdaptiveDT(Real & dt, bool allowToGrow = true, bool allowToShrink = true) override;

  Real localErrorEstimate();

  const Real _tolerance;
  const Real _Kint;
  const Real _Kpro;
  const Real _Kder;
  const int _start_adapting;
  const int _local_error_order;
  NonlinearSystemBase & _nl;
  AuxiliarySystem & _aux;
  const bool _check_aux_error;

  Real _cutback_factor;

  Real _local_error_estimate_current;
  Real _local_error_estimate_old;
  Real _local_error_estimate_older;
  Real _local_error_estimate_older_save;
  Real _R;
  
  std::vector<std::string> _check_variables;

// private:
//   Real localErrorEstimate1(const NumericVector<Number> & s1, NumericVector<Number> & s2);
//   Real localErrorEstimate2(const NumericVector<Number> & s1, NumericVector<Number> & s2, NumericVector<Number> & s3);
};
