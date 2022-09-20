[Mesh]
  type = GeneratedMesh
  dim = 1
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [dt]
    type = TimeDerivative
    variable = u
  []
[]

[Executioner]
  type = Transient
  end_time = 20.0
  verbose = true

  [TimeStepper]
    type = PIDAdaptiveDT
    dt = 1.0
    tolerance = 0.001
    K_integral = 0.05
    K_proportional = 0.05
    K_derivative = 0.05
  []
[]

[Postprocessors]
  [_dt]
    type = TimestepSize
  []
[]

[Outputs]
  csv = true
[]
