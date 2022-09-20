//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TimestepSize.h"
#include "FEProblem.h"

registerMooseObject("MooseApp", TimestepSize);

InputParameters
TimestepSize::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription("Reports the timestep size");
  params.addParam<bool>("next_timestep", false, "Get the next (unconstrained) time step predicted by time stepper");
  return params;
}

TimestepSize::TimestepSize(const InputParameters & parameters)
  : GeneralPostprocessor(parameters), _feproblem(dynamic_cast<FEProblemBase &>(_subproblem)),
    _isnextTstepper(getParam<bool>("next_timestep"))
{
  if (_isnextTstepper)
    _timestepper = dynamic_cast<Transient *>(_app.getExecutioner())->getTimeStepper();
}

Real
TimestepSize::getValue()
{
  if (_isnextTstepper)
    return _timestepper->getNextDT();
  else
    return _feproblem.dt();
}
