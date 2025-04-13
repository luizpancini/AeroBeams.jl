## Index

```@index
Pages = ["publicAPI.md"]
```

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
update_Airfoil_params!
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
remove_all_springs_from_beams!
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

### Creating a hinge axis constraint
```@docs
create_HingeAxisConstraint
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

### Creating a sample model
```@docs
geometrical_properties_Pazy
tip_loss_function_Pazy
typical_section_data
create_Pazy
create_PazyFFWT
create_SMW
create_TDWing
create_Helios
create_conventional_HALE
create_BWB
create_HealyFFWT
create_HealyBaselineFFWT
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
axial
isotropic_stiffness_matrix
inertia_matrix
rotation_tensor_E321
rotation_tensor_E213
rotation_tensor_E231
rotation_tensor_E313
rotation_tensor_WM
scaled_rotation_parameters
rotation_angle
rotation_angle_limited
rotation_parameters_WM
rotation_parameters_Rodrigues
yrp_from_rotation_tensor
quaternion_from_rotation_tensor
WM_to_yrp
yrp_to_WM
rotation_between_WM
mode_tracking
get_FFT_and_PSD
Newton_solver
backward_extrapolation
NASA_frames_loader
GU_frames_loader
```

### Visualizing the results
```@docs
plot_undeformed_assembly
plot_steady_deformation
plot_steady_outputs
plot_mode_shapes
plot_dynamic_deformation
plot_time_outputs
plot_snapshots
```