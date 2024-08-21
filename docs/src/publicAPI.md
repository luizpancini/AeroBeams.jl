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
create_Airfoil
create_flapped_Airfoil
```

### Creating an aerodynamic surface
```@docs
create_AeroSurface
```

### Creating an atmosphere
```@docs
standard_atmosphere
```

### Creating boundary conditions
```@docs
create_BC
```

### Creating a beam
```@docs
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

### Creating and solving a problem
```@docs
InitialVelocitiesUpdateOptions
create_SteadyProblem
create_TrimProblem
create_EigenProblem
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
tip_loss_factor_Pazy
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
isotropic_stiffness_matrix
inertia_matrix
rotation_tensor_E321
rotation_tensor_E313
rotation_tensor_WM
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