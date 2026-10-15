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
create_FlapParameters
create_Airfoil
update_all_airfoil_parameters!
modified_database
```

The following airfoils are predefined: `FlatPlate`, `NACA0002`, `NACA0012`, `NACA0015`, `NACA0018`, `NACA23012A`, `VERTOL23010`, `HeliosWingAirfoil`, `HeliosPodAirfoil`, `BWBAirfoil` and `cHALEairfoil`. Use `deepcopy` to get an independent copy before modifying one of them.

### Creating an aerodynamic surface
```@docs
create_AeroSurface
```

### Creating an atmosphere
```@docs
standard_atmosphere
air_properties_from_pressure_and_temperature
```

### Creating boundary conditions
```@docs
create_BC
reset_BC_beam!
```

### Creating a beam
```@docs
create_Beam
update_beam_orientation!
add_distributed_loads_to_beam!
add_initial_conditions_to_beam!
add_point_inertia_to_beam!
add_spring_to_beam!
add_spring_to_beams!
add_surface_to_beam!
remove_all_loads_from_beams!
remove_all_initial_conditions_from_beams!
remove_all_springs_from_beams!
remove_all_point_inertias_from_beams!
remove_surface_from_beams!
```

### Creating a gust
```@docs
create_SharpEdgedGust
create_OneMinusCosineGust
create_Continuous1DGust
create_DiscreteSpaceGust
create_Continuous1DSpaceGust
create_Continuous2DSpaceGust
stochastic_gust_velocity_from_white_noise
remove_gust!
```

### Creating a link
```@docs
create_TrimLoadsLink
create_FlapLink
```

### Creating a model
```@docs
create_Model
add_BCs_to_model!
update_model_BCs!
update_model_links!
update_model_hinge_axis_constraints!
update_initial_rotation_of_basis_A!
update_u_A!
update_v_A!
update_vdot_A!
update_ω_A!
update_ωdot_A!
```

### Creating a point inertia
```@docs
create_PointInertia
```

### Creating a hinge axis constraint
```@docs
create_HingeAxisConstraint
```

### Creating and solving a problem
```@docs
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
create_SMW
create_TDWing
create_Helios
create_cHALE
create_BWB
create_HealyBaselineFFWT
create_HealySideslipFFWT
create_HealyLCOFFWT
tip_loss_function_HealyLCOFFWT
```

### Creating a spring
```@docs
create_Spring
```

### Creating a system solver
```@docs
create_LinearSolver
create_NewtonRaphson
create_NLsolver
create_ArclengthSolver
create_TrimSolver
```

### Creating a units system
```@docs
create_UnitsSystem
feet2meters
pounds2kilograms
rads2Hz
Hz2rads
rpm2rads
rads2rpm
```

### Utilities
```@docs
round_off
round_off!
AeroBeams.average
rms
tilde
axial
gauss_legendre
isotropic_stiffness_matrix
anisotropic_stiffness_matrix
inertia_matrix
rotation_tensor_E321
rotation_tensor_E213
rotation_tensor_E231
rotation_tensor_E313
rotation_tensor_quaternion
transform_to_WM
rotation_tensor_WM
scaled_rotation_parameters
twist_angle
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
filter_tracked_modes
find_flutter
find_divergence
blockdiag_dense
moving_average
get_FFT_and_PSD
local_extrema
backward_extrapolation
count_of_exceedance
frequency_of_exceedance
gcd_floats_approx
spectrogram
NASA_frames_loader
GU_frames_loader
```

### Visualizing the results
```@docs
get_undeformed_airfoil_coords
plot_undeformed_assembly
plot_steady_deformation
plot_steady_deformations
plot_steady_outputs
plot_mode_shapes
plot_mode_shapes_animation
plot_dynamic_deformation
plot_dynamic_deformations
plot_time_outputs
plot_snapshots
```