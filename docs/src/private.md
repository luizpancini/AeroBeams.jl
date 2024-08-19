## Private API

### Aerodynamics
```@docs
AeroBeams.aero_steady_kinematics!
AeroBeams.aero_unsteady_kinematics!
AeroBeams.nondimensional_flow_parameters!
AeroBeams.local_gust_velocity!
AeroBeams.flap_deflection_rates!
AeroBeams.aero_coefficients!
AeroBeams.aero_state_matrices!
AeroBeams.aero_loads_resultants!
AeroBeams.attached_flow_aero_coefficients!
AeroBeams.update_airfoil_parameters!
AeroBeams.effective_angle_of_attack!
AeroBeams.attached_flow_cn!
AeroBeams.attached_flow_cm!
AeroBeams.attached_flow_ct!
AeroBeams.pitch_plunge_effective_normalwash
AeroBeams.flap_effective_normalwash
AeroBeams.flap_normalwash
AeroBeams.gust_effective_normalwash
AeroBeams.flap_normalwash_rate
AeroBeams.cnαUₙTQC_rate
AeroBeams.attached_flow_state_matrices!
AeroBeams.BLi_aero_coefficients!
AeroBeams.BLi_kinematics!
AeroBeams.BLi_nonlinear_states!
AeroBeams.BLi_motion_qualifiers!
AeroBeams.BLi_breakpoint_angles!
AeroBeams.BLi_time_delays!
AeroBeams.BLi_separation_points!
AeroBeams.BLi_quasi_steady_separation_points
AeroBeams.BLi_stall_time!
AeroBeams.BLi_DSV_loads!
AeroBeams.BLi_cn!
AeroBeams.BLi_cm!
AeroBeams.BLi_ct!
AeroBeams.BLi_state_matrices!
AeroBeams.BLo_aero_coefficients!
AeroBeams.BLo_nonlinear_states!
AeroBeams.BLo_motion_qualifiers!
AeroBeams.BLo_breakpoint_angle!
AeroBeams.BLo_separation_points!
AeroBeams.BLo_stall_time!
AeroBeams.BLo_time_delays!
AeroBeams.BLo_update_impulsive_parameters!
AeroBeams.BLo_vortex_accumulation_rate!
AeroBeams.BLo_cn!
AeroBeams.BLo_cm!
AeroBeams.BLo_ct!
AeroBeams.BLo_state_matrices!
AeroBeams.update_initial_aero_states!
```

### AeroProperties
```@docs
AeroBeams.FlowParameters
AeroBeams.FlowAnglesAndRates
AeroBeams.FlowVelocitiesAndRates
AeroBeams.AeroCoefficients
AeroBeams.BLiNamedStates
AeroBeams.BLoNamedStates
AeroBeams.BLiKinematics
AeroBeams.BLiFlowVariables
AeroBeams.BLoFlowVariables
AeroBeams.BLiComplementaryVariables
AeroBeams.BLoComplementaryVariables
AeroBeams.AeroVariables
AeroBeams.AeroProperties
AeroBeams.initial_F_χ_χ
AeroBeams.initial_F_χ_Vdot
AeroBeams.initial_F_χ_Ωdot
```

### AeroSolver
```@docs
AeroBeams.theodorsen_flap_constants
AeroBeams.inflow_arrays
```

### Airfoil
```@docs
AeroBeams.AttachedFlowParameters
AeroBeams.BLiParameters
AeroBeams.BLoParameters
AeroBeams.get_airfoil_coordinates
```

### Atmosphere
```@docs
AeroBeams.air_properties_from_pressure_and_temperature
```

### Boundary conditions
```@docs
AeroBeams.update_BC_data!
```

### Beam
```@docs
AeroBeams.validate_beam!
AeroBeams.validate_sectional_matrices
AeroBeams.validate_rotation_parametrization
AeroBeams.validate_connected_beams
AeroBeams.validate_initial_conditions!
AeroBeams.validate_normalized_nodal_positions!
AeroBeams.validate_hinged_nodes!
AeroBeams.validate_distributed_loads!
AeroBeams.velocity_dofs_to_update!
AeroBeams.get_rotation_tensor!
AeroBeams.initial_displacements_derivatives!
AeroBeams.create_beam_elements!
AeroBeams.set_nodal_coordinates!
```

### Core
```@docs
AeroBeams.element_arrays!
AeroBeams.special_node_arrays!
AeroBeams.element_velocities_basis_b!
AeroBeams.element_accelerations_basis_b!
AeroBeams.element_states!
AeroBeams.element_states_rates!
AeroBeams.element_rotation_variables!
AeroBeams.element_distributed_loads!
AeroBeams.gravitational_loads!
AeroBeams.distributed_external_loads!
AeroBeams.aerodynamic_loads!
AeroBeams.wrapper_aerodynamic_loads_from_states!
AeroBeams.wrapper_aerodynamic_loads_from_states_rates!
AeroBeams.aero_loads_core!
AeroBeams.interpolate_distributed_loads
AeroBeams.element_strains!
AeroBeams.element_momenta!
AeroBeams.element_momenta_rates!
AeroBeams.element_residual!
AeroBeams.distributed_loads_derivatives_rotation_parameters!
AeroBeams.gravitational_loads_derivatives_rotation_parameters
AeroBeams.distributed_external_loads_derivatives_rotation_parameters
AeroBeams.aero_loads_derivatives_rotation_parameters
AeroBeams.aero_derivatives!
AeroBeams.element_jacobian!
AeroBeams.element_inertia!
AeroBeams.element_nodal_states!
AeroBeams.special_node_states!
AeroBeams.spring_loads!
AeroBeams.special_node_residual!
AeroBeams.special_node_follower_loads_derivatives_rotation_parameters!
AeroBeams.special_node_jacobian!
AeroBeams.update_special_node_jacobian!
AeroBeams.spring_loads_jacobians!
AeroBeams.element_modal_states
AeroBeams.update_states!
AeroBeams.reset_dual_numbers
AeroBeams.convert_to_values
```

### Element
```@docs
AeroBeams.ElementalStates
AeroBeams.ComplementaryElementalStates
AeroBeams.ElementalStatesRates
AeroBeams.ComplementaryElementalStatesRates
AeroBeams.NodalStates
AeroBeams.Element
AeroBeams.get_hinged_nodes_matrices
AeroBeams.get_element_distributed_loads
AeroBeams.update_element_distributed_loads!
AeroBeams.add_point_inertia_to_element!
```

### Gust
```@docs
AeroBeams.SharpEdgedGust
AeroBeams.OneMinusCosineGust
AeroBeams.Continuous1DGust
AeroBeams.DiscreteSpaceGust
AeroBeams.Continuous1DSpaceGust
AeroBeams.Continuous2DSpaceGust
AeroBeams.stochastic_gust_velocity_from_white_noise
```

### Links
```@docs
AeroBeams.TrimLoadsLink
AeroBeams.FlapLink
```

### Model
```@docs
AeroBeams.Model
AeroBeams.validate_model!
AeroBeams.validate_and_update_motion_basis_A!
AeroBeams.assemble_model!
AeroBeams.inertia_properties!
AeroBeams.set_atmosphere!
AeroBeams.update_linked_flap_deflections!
AeroBeams.update_number_gust_states!
AeroBeams.update_loads_trim_links_global_ids!
AeroBeams.update_spring_nodes_ids!
AeroBeams.update_relative_rotation_constraint_elements_ids!
AeroBeams.initialize_basis_A_rotation!
AeroBeams.update_initial_conditions!
AeroBeams.set_BCs!
AeroBeams.get_special_nodes!
AeroBeams.get_system_indices!
AeroBeams.update_relative_rotation_constraint_data!
```

### Problem
```@docs
AeroBeams.Problem
AeroBeams.ModeShape
AeroBeams.set_initial_states!
AeroBeams.initialize_system_arrays!
AeroBeams.precompute_distributed_loads!
AeroBeams.solve_steady!
AeroBeams.get_mode_shapes!
AeroBeams.fix_nodal_modal_states!
AeroBeams.solve_dynamic!
AeroBeams.initialize_time_variables!
AeroBeams.solve_initial_dynamic!
AeroBeams.copy_initial_states!
AeroBeams.update_initial_velocities!
AeroBeams.time_march!
AeroBeams.adaptable_time_march!
AeroBeams.copy_state
AeroBeams.restore_state!
AeroBeams.update_time_variables!
AeroBeams.update_basis_A_orientation!
AeroBeams.get_equivalent_states_rates!
AeroBeams.update_BL_complementary_variables!
AeroBeams.BL_stall_boundary
AeroBeams.solve_time_step!
AeroBeams.save_time_step_data!
```

### Rotation constraint
```@docs
AeroBeams.RotationConstraint
```

### Sample models
```@docs
AeroBeams.Pazy_nodal_positions
AeroBeams.Pazy_stiffness_matrices
AeroBeams.Pazy_inertia_matrices
```

### Special node
```@docs
AeroBeams.SpecialNode
```

### Spring
```@docs
AeroBeams.Spring
```

### System solver
```@docs
AeroBeams.SystemSolver
AeroBeams.NewtonRaphson
AeroBeams.solve_NewtonRaphson!
AeroBeams.assemble_system_arrays!
AeroBeams.solve_linear_system!
AeroBeams.line_search
AeroBeams.line_search_step_size
AeroBeams.save_load_factor_data!
AeroBeams.update_maximum_aero_loads!
```

### Units system
```@docs
AeroBeams.UnitsSystem
AeroBeams.validate_units_system
```

### Utilities
```@docs
AeroBeams.a1
AeroBeams.a2
AeroBeams.a3
AeroBeams.I3
AeroBeams.I6
AeroBeams.divide_inplace!
AeroBeams.multiply_inplace!
AeroBeams.curvature_quantities
AeroBeams.position_vector_from_curvature
AeroBeams.rotation_tensor_from_curvature
AeroBeams.rotation_parameter_scaling
AeroBeams.rotation_tensor_derivatives_scaled_parameters
AeroBeams.scaling_derivatives_extended_parameters
AeroBeams.rotation_tensor_derivatives_extended_parameters
AeroBeams.rotation_tensor_time_derivative
AeroBeams.rotation_tensor_derivatives_time_extended_parameters
AeroBeams.tangent_operator_transpose_inverse_WM
AeroBeams.tangent_tensor_functions_derivatives_extended_parameters
AeroBeams.force_scaling
AeroBeams.rotation_angle
AeroBeams.highest_in_rowcol
```

### Visualization
```@docs
AeroBeams.plot_output_of_x1
AeroBeams.plot_output_of_time!
AeroBeams.get_undeformed_airfoil_coords
AeroBeams.plot_BCs!
AeroBeams.draw_BC!
AeroBeams.draw_concentrated_force!
AeroBeams.draw_concentrated_moment!
AeroBeams.draw_generalized_displacement!
AeroBeams.plot_distributed_loads!
AeroBeams.draw_distributed_forces!
AeroBeams.draw_distributed_moments!
AeroBeams.draw_aero_loads!
AeroBeams.draw_circular_vector!
```