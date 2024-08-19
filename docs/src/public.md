## Public API

### Creating an aerodynamic derivatives solver
```@docs
AD
FD
```

### Creating an aerodynamic solver
```@docs
QuasiSteady
Indicial
BLi
BLo
Inflow
TableLookup
ThinAirfoilTheory
IndicialGust
```

### Creating an airfoil
```@docs
Airfoil
create_Airfoil
create_flapped_Airfoil
```

### Creating an aerodynamic surface
```@docs
AeroSurface
create_AeroSurface
```

### Creating an atmosphere
```@docs
Atmosphere
standard_atmosphere
```

### Creating boundary conditions
```@docs
BC
create_BC
```

### Creating a beam
```@docs
Beam
create_Beam
update_beam!
add_point_inertias_to_beam!
add_loads_to_beam!
add_initial_displacements_and_velocities_to_beam!
add_springs_to_beam!
add_spring_to_beams!
```

### Creating a gust
```@docs
create_SharpEdgedGust
create_OneMinusCosineGust
create_Continuous1DGust
create_DiscreteSpaceGust
create_Continuous1DSpaceGust
create_Continuous2DSpaceGust
```

### Creating a link
```@docs
create_TrimLoadsLink
create_FlapLink
```

### Creating a model
```@docs
create_Model
update_model!
set_motion_basis_A!
```

### Creating a point inertia
```@docs
PointInertia
```

### Creating a problem
```@docs
InitialVelocitiesUpdateOptions
SteadyProblem
create_SteadyProblem
TrimProblem
create_TrimProblem
EigenProblem
create_EigenProblem
DynamicProblem
create_DynamicProblem
solve!
solve_eigen!
```

### Creating a rotation constraint
```@docs
create_RotationConstraint
```

### Creating a sample model
```@docs
create_Pazy
create_PazyFFWT
Pazy_tip_loss_factor
create_Helios
create_conventional_HALE
create_BWB
```

### Creating a spring
```@docs
create_Spring
```

### Creating a system solver
```@docs
create_NewtonRaphson
```

### Creating a units system
```@docs
create_UnitsSystem
```

### Utilities
```@docs
round_off!
rms
tilde
mul3
isotropic_stiffness_matrix
inertia_matrix
rotation_tensor_E321
rotation_tensor_E313
rotation_tensor_WM
tangent_operator_transpose_WM
tangent_tensor_transpose_derivatives_extended_parameters
rotation_parameters_WM
ypr_from_rotation_tensor
quaternion_from_rotation_tensor
mode_tracking
get_FFT_and_PSD
```

### Visualizing the results
```@docs
plot_undeformed_assembly
plot_steady_deformation
plot_steady_outputs
plot_mode_shapes
plot_dynamic_deformation
plot_time_outputs
```