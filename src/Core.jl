
show_debug() = true

# Calls f on every element, in order, looping over one group of concretely-typed elements at a time
@inline foreach_element(f::F, ::Tuple{}) where {F} = nothing
@inline function foreach_element(f::F, elementGroups::Tuple) where {F}
    for element in first(elementGroups)
        f(element)
    end
    foreach_element(f, Base.tail(elementGroups))
end

# Assembles the residual vector of the system of equations
function f!(problem)

    @unpack x,residual = problem.systemArrays
    @unpack σ = problem.systemSolver.loadStepState
    @unpack model = problem
    @unpack beams,hingeAxisConstraints = model
    @unpack elementGroups,specialNodes = model.assemblyData
    timeIndex = time_index(problem)
    timeNow = time_now(problem)
    Δt = time_step(problem)
    isDynamicProblem = problem isa DynamicProblem

    residual!(residual,model,beams,elementGroups,specialNodes,hingeAxisConstraints,x,σ,timeIndex,timeNow,Δt,isDynamicProblem)

    return residual

end

# Assembles the Jacobian matrix of the system of equations
function j!(problem)

    @unpack x,jacobian = problem.systemArrays
    @unpack σ = problem.systemSolver.loadStepState
    @unpack model = problem
    @unpack beams,hingeAxisConstraints = model
    @unpack elementGroups,specialNodes = model.assemblyData
    skipAeroJacobianUpdate = get_skipAeroJacobianUpdate(problem.systemSolver)
    timeNow = time_now(problem)
    Δt = time_step(problem)
    isDynamicProblem = problem isa DynamicProblem

    jacobian!(jacobian,model,elementGroups,specialNodes,hingeAxisConstraints,x,σ,skipAeroJacobianUpdate,timeNow,Δt,isDynamicProblem)

    return jacobian

end

# Assembles the inertia matrix (only for eigenproblems and trim problems, when required)
function i!(problem)

    @unpack inertia = problem.systemArrays
    @unpack model = problem
    @unpack elementGroups = model.assemblyData
    timeNow = time_now(problem)
    Δt = time_step(problem)
    isDynamicProblem = problem isa DynamicProblem

    # Loop over elements and set their contributions to the inertia matrix
    foreach_element(elementGroups) do element
        element_inertia!(inertia,model,element,timeNow,Δt,isDynamicProblem)
    end

end
                        
# Assembles the residual vector and Jacobian matrix of the system of equations
function fj!(problem)

    @unpack x,residual,jacobian = problem.systemArrays
    @unpack σ = problem.systemSolver.loadStepState
    @unpack model = problem
    @unpack beams,hingeAxisConstraints = model
    @unpack elementGroups,specialNodes = model.assemblyData
    skipJacobianUpdate = get_skipJacobianUpdate(problem.systemSolver)
    skipAeroJacobianUpdate = get_skipAeroJacobianUpdate(problem.systemSolver)
    timeIndex = time_index(problem)
    timeNow = time_now(problem)
    Δt = time_step(problem)
    isDynamicProblem = problem isa DynamicProblem

    # Set residual vector
    residual!(residual,model,beams,elementGroups,specialNodes,hingeAxisConstraints,x,σ,timeIndex,timeNow,Δt,isDynamicProblem)

    # Set Jacobian matrix, if applicable
    skipJacobianUpdate || jacobian!(jacobian,model,elementGroups,specialNodes,hingeAxisConstraints,x,σ,skipAeroJacobianUpdate,timeNow,Δt,isDynamicProblem)

    return residual,jacobian

end

# Updates the residual vector of the system of equations
function residual!(residual,model,beams,elementGroups,specialNodes,hingeAxisConstraints,x,σ,timeIndex,timeNow,Δt,isDynamicProblem)

    # Set flap deflections and their rates
    for beam in beams
        flap_deflection!(beam,x,timeIndex)
    end

    # Set contributions from the elements
    foreach_element(elementGroups) do element
        element_residual!(residual,model,element,x,σ,timeIndex,timeNow,Δt,isDynamicProblem)
    end

    # Update states of all special nodes first (necessary for convergence with doubly-attached springs)
    for specialNode in specialNodes
        special_node_states!(model,specialNode,x,σ)
    end

    # Set contributions from the special nodes
    for specialNode in specialNodes
        special_node_residual!(residual,model,specialNode)
    end

    # Set contributions from the hinge axis constraints
    for constraint in hingeAxisConstraints
        hac_residual!(residual,model,constraint,x)
    end

end

# Updates the Jacobian matrix of the system of equations
function jacobian!(jacobian,model,elementGroups,specialNodes,hingeAxisConstraints,x,σ,skipAeroJacobianUpdate,timeNow,Δt,isDynamicProblem)

    # Reset Jacobians associated with spring loads and hinge axis constraints, and the trim flap deflections' columns
    reset_jacobians!(jacobian,model,specialNodes,hingeAxisConstraints)

    # Set contributions from the elements
    foreach_element(elementGroups) do element
        element_jacobian!(jacobian,model,element,timeNow,Δt,isDynamicProblem,skipAeroJacobianUpdate)
    end

    # Set contributions from the special nodes
    for specialNode in specialNodes
        special_node_jacobian!(jacobian,model,specialNode,x,σ)
    end

    # Set contributions from the hinge axis constraints
    for constraint in hingeAxisConstraints
        hac_jacobian!(jacobian,model,constraint)
    end

end

# Updates the flap deflections of beam surfaces and their rates
function flap_deflection!(beam,x,timeIndex)

    # Unpack surface
    @unpack surface = beam.attachments

    # Skip beams without aerodynamic surfaces
    isnothing(surface) && return

    # Current flap deflection (either a trim state or a time-dependent value) and its rates
    if surface.assemblyData.δIsTrimVariable
        @unpack DOF_δ,δMultiplier = surface.assemblyData
        δ = x[DOF_δ[1]] * δMultiplier
        δdot = δddot = 0.0
    else
        @unpack δVec,δdotVec,δddotVec = surface.assemblyData
        δ = δVec[timeIndex]
        δdot = δdotVec[timeIndex]
        δddot = δddotVec[timeIndex]
    end
    
    # Update surface state
    setfield!(surface.state, :δ, δ)
    setfield!(surface.state, :δdot, δdot)
    setfield!(surface.state, :δddot, δddot)

end

# Computes the contributions from the current element to the residual array
function element_residual!(residual,model,element,x,σ,timeIndex,timeNow,Δt,isDynamicProblem)

    @unpack invForceScaling = model.assemblyData
    @unpack Δℓ,R0,R0T,hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat = element.geometry
    @unpack eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,eqs_Fχ,isSpecialNode1,isSpecialNode2,eqsNode1Set,eqsNode2Set = element.assemblyIndices
    @unpack u,p,F,M,V,Ω,γ,κ,P,H,χ = element.state.midpoint
    @unpack udot,pdot,Pdot,Hdot,χdot = element.state.rates
    @unpack v,ω = element.stateCache.kinematicsBasisb
    @unpack RR0,RR0T,RdotR0,HT,HTinv = element.stateCache.rotationMatrices
    @unpack f1,f2,m1,m2 = element.stateCache.loadsResultants
    @unpack A,B = element.aero.stateCache.stateMatrices
    @unpack F_u1,F_u2,F_p1,F_p2,F_F1,F_F2,F_M1,F_M2,F_V,F_Ω,F_χ = element.stateCache.residuals

    # Generalized velocities of basis b at element's midpoint, resolved in basis A
    element_velocities_basis_b!(model,element,σ)

    # States
    element_states!(model,element,x)

    # States' rates
    isDynamicProblem && element_states_rates!(element,Δt)
    
    # Rotation variables
    element_rotation_variables!(element,isDynamicProblem)
    
    # Distributed loads
    element_distributed_loads!(model,element,σ,timeIndex,timeNow,Δt,isDynamicProblem)

    # Static terms
    element_static_residual!(F_u1,F_u2,F_p1,F_p2,F_F1,F_F2,F_M1,F_M2,Δℓ,R0,u,p,F,M,γ,κ,RR0,HTinv,f1,f2,m1,m2)

    # Steady terms
    element_steady_residual!(F_u1,F_u2,F_p1,F_p2,F_V,F_Ω,F_χ,Δℓ,u,V,Ω,P,H,χ,v,ω,RR0,RR0T,A,B)

    # Transient dynamic terms
    isDynamicProblem && element_transient_residual!(F_u1,F_u2,F_p1,F_p2,F_V,F_Ω,F_χ,Δℓ,R0T,P,H,udot,pdot,Pdot,Hdot,χdot,RR0,RdotR0,HT)

    # Insert resultants into the residual array
    insert_element_residual!(residual,F_u1,F_u2,F_p1,F_p2,F_F1,F_F2,F_M1,F_M2,F_V,F_Ω,F_χ,invForceScaling,eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,eqs_Fχ,isSpecialNode1,isSpecialNode2,eqsNode1Set,eqsNode2Set,hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat)

    # @logmsg LogLevel(show_debug() ? 0 : -1) "ELEMENT $(element.assemblyData.globalID) RESIDUALS"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_u1: $F_u1"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p1: $F_p1"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_F1: $F_F1"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_M1: $F_M1"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_V: $F_V"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_Ω: $F_Ω"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_χ: $F_χ"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "norm(residual): $(norm(residual))"

    return nothing
    
end

# Computes the generalized velocities of basis b at the element's midpoint, resolved in basis A
function element_velocities_basis_b!(model,element,σ)

    @unpack v_A,ω_A,R_A,R_AT = model.kinematicsBasisA.state
    @unpack r = element.assemblyData
    @unpack v,ω = element.stateCache.kinematicsBasisb

    # Translational velocity
    v .= σ * SVector(R_AT * (v_A + cross(ω_A, R_A*r)))

    # Angular velocity
    ω .= σ * SVector(R_AT * ω_A)

end

# Computes the generalized accelerations of basis b at the element's midpoint, resolved in basis A
function element_accelerations_basis_b!(model,element,σ)

    @unpack v_A,ω_A,vdot_A,ωdot_A,R_A,R_AT = model.kinematicsBasisA.state
    @unpack r = element.assemblyData
    @unpack vdot,ωdot = element.stateCache.kinematicsBasisb

    # Translational acceleration
    vdot .= σ * SVector(R_AT * (vdot_A + 
                        cross(ωdot_A, R_A*r) + 
                        cross(ω_A, cross(ω_A, R_A*r)) - 
                        cross(ω_A, R_A*v_A)
                        ))

    # Angular acceleration
    ωdot .= σ * SVector(R_AT * ωdot_A)

end

# Sets the midpoint states (generalized displacements, forces, strains, velocities and aerodynamic) of the element
function element_states!(model,element,x)

    @unpack forceScaling = model.assemblyData
    @unpack DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
    @unpack C11,C12,C21,C22,I11,I12,I21,I22 = element.sectional
    @unpack u,p,F,M,V,Ω,γ,κ,P,H,χ = element.state.midpoint

    # Actual midpoint states
    u .= x[DOF_u]
    p .= x[DOF_p]
    F .= x[DOF_F] * forceScaling
    M .= x[DOF_M] * forceScaling
    V .= x[DOF_V]
    Ω .= x[DOF_Ω]
    χ .= x[DOF_χ]

    # Complementary midpoint states
    γ .= C11 * F + C12 * M
    κ .= C21 * F + C22 * M
    P .= I11 * V + I12 * Ω
    H .= I21 * V + I22 * Ω

    # @logmsg LogLevel(show_debug() ? 0 : -1) "ELEMENT $(element.assemblyData.globalID) STATES"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "u: $u"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "p: $p"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F: $F"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "M: $M"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "V: $V"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "Ω: $Ω"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "χ: $χ"

end

# Sets the states' rates of the element
function element_states_rates!(element,Δt)

    @unpack I11,I12,I21,I22 = element.sectional
    @unpack u,p,V,Ω,χ = element.state.midpoint
    @unpack udot,pdot,Vdot,Ωdot,Pdot,Hdot,χdot = element.state.rates
    @unpack udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv,χdotEquiv = element.stateCache.equivStatesRates
    
    # Set rates for finite time step
    if isfinite(Δt)
        udot .= 2/Δt * u - udotEquiv
        pdot .= 2/Δt * p - pdotEquiv
        Vdot .= 2/Δt * V - VdotEquiv
        Ωdot .= 2/Δt * Ω - ΩdotEquiv
        χdot .= 2/Δt * χ - χdotEquiv
        Pdot .= I11 * Vdot + I12 * Ωdot
        Hdot .= I21 * Vdot + I22 * Ωdot
    end

end

# Computes the rotation variables for the element (values only, derivatives are computed in element_rotation_derivatives!)
function element_rotation_variables!(element,isDynamicProblem)

    @unpack R0 = element.geometry
    @unpack p = element.state.midpoint
    @unpack pdot = element.state.rates
    @unpack R,RR0,RR0T,RdotR0,HT,HTinv = element.stateCache.rotationMatrices

    # Set rotation tensor from basis b to basis B, resolved in basis A, and get associated variables
    ps,ps0,υ²,pspsT,ps_p = R_WM!(R,p)
    
    # Rotation tensor from basis A to basis B, resolved in basis A, and its transpose
    RR0 .= R * R0 
    RR0T .= RR0'
    
    # Functions of tangent operator tensor (HT and its inverse)
    HT_WM!(HT,ps,ps0,υ²,pspsT,ps_p)
    HTinv_WM!(HTinv,ps,ps0,pspsT,ps_p)
    
    # Time derivative of rotation tensor
    isDynamicProblem && RdotR0!(RdotR0,R,R0,HT,pdot)
    
end

# Computes the nodal resultants from distributed loads on the current element, resolved in basis A
function element_distributed_loads!(model,element,σ,timeIndex,timeNow,Δt,isDynamicProblem)

    @unpack f1,f2,m1,m2 = element.stateCache.loadsResultants
    @unpack f1_ext,f2_ext,m1_ext,m2_ext,ff1_A,ff2_A,mf1_A,mf2_A,ff1_b,ff2_b,mf1_b,mf2_b,f_g,m_g = element.stateCache.structuralLoadsResultants
    @unpack f1χ,f2χ,m1χ,m2χ = element.aero.stateCache.loadsResultants
    @unpack Ac,Ai,An,Af,Ag,A,B = element.aero.stateCache.stateMatrices

    # Gravitational loads 
    gravitational_loads!(f_g,m_g,model,element,σ)

    # Externally applied distributed loads
    distributed_external_loads!(f1_ext,f2_ext,m1_ext,m2_ext,ff1_A,ff2_A,mf1_A,mf2_A,ff1_b,ff2_b,mf1_b,mf2_b,element,σ,timeIndex)

    # Aerodynamic loads (and state arrays)
    if !isnothing(element.aero.surface)
        aerodynamic_loads_and_arrays!(f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,A,B,model,element,timeNow,Δt,isDynamicProblem)
    end

    # Total nodal resultants from distributed loads
    f1 .= f_g + f1_ext + f1χ
    f2 .= f_g + f2_ext + f2χ
    m1 .= m_g + m1_ext + m1χ
    m2 .= m_g + m2_ext + m2χ

end

# Computes the nodal resultants from the distributed gravitational loads on the current element
function gravitational_loads!(f_g,m_g,model,element,σ)

    @unpack gravityVector = model
    @unpack R_AT = model.kinematicsBasisA.state
    @unpack Δℓ = element.geometry
    @unpack ηtilde,μ = element.sectional
    @unpack RR0,RR0T = element.stateCache.rotationMatrices

    # Nodal distributed weight force and moment vectors, resolved in basis A
    if any(!iszero(gravityVector)) 
        f_g .= σ * Δℓ/2 * μ * R_AT * gravityVector
        m_g .= SMatrix(RR0 * ηtilde) * SVector(RR0T * f_g)
    end

end

# Computes the nodal resultants from the externally applied distributed loads on the current element, all resolved in basis A
function distributed_external_loads!(f1_ext,f2_ext,m1_ext,m2_ext,ff1_A,ff2_A,mf1_A,mf2_A,ff1_b,ff2_b,mf1_b,mf2_b,element,σ,timeIndex)

    @unpack f_A,m_A,f_b,m_b,ff_A,mf_A,ff_b,mf_b,has_f_A,has_m_A,has_f_b,has_m_b,has_ff_A,has_mf_A,has_ff_b,has_mf_b = element.distLoads
    @unpack R0 = element.geometry
    @unpack R,RR0 = element.stateCache.rotationMatrices

    # Reset to zero
    zero_arrays!(f1_ext,f2_ext,m1_ext,m2_ext) 

    ## Dead loads
    # --------------------------------------------------------------------------
    # Add dead forces initially resolved in basis A
    if has_f_A
        f1_ext .+= σ * f_A[1,:,timeIndex]
        f2_ext .+= σ * f_A[2,:,timeIndex]
    end

    # Add dead moments initially resolved in basis A
    if has_m_A
        m1_ext .+= σ * m_A[1,:,timeIndex]
        m2_ext .+= σ * m_A[2,:,timeIndex]
    end

    # Add dead forces initially resolved in basis b
    if has_f_b
        f1_ext .+= R0 * σ * f_b[1,:,timeIndex]
        f2_ext .+= R0 * σ * f_b[2,:,timeIndex]
    end

    # Add dead moments initially resolved in basis b
    if has_m_b
        m1_ext .+= R0 * σ * m_b[1,:,timeIndex]
        m2_ext .+= R0 * σ * m_b[2,:,timeIndex]
    end

    ## Follower loads
    # --------------------------------------------------------------------------
    # Set and add follower forces initially resolved in basis A
    if has_ff_A
        # Set values
        @views ff1_A .= σ * ff_A[1,:,timeIndex]
        @views ff2_A .= σ * ff_A[2,:,timeIndex]
        # Add
        f1_ext .+= R * ff1_A
        f2_ext .+= R * ff2_A
    end

    # Set and add follower moments initially resolved in basis A
    if has_mf_A
        # Set values
        @views mf1_A .= σ * mf_A[1,:,timeIndex]
        @views mf2_A .= σ * mf_A[2,:,timeIndex]
        # Add
        m1_ext .+= R * mf1_A
        m2_ext .+= R * mf2_A
    end

    # Set and add follower forces initially resolved in basis b
    if has_ff_b
        # Set values
        @views ff1_b .= σ * ff_b[1,:,timeIndex]
        @views ff2_b .= σ * ff_b[2,:,timeIndex]
        # Add
        f1_ext .+= RR0 * ff1_b
        f2_ext .+= RR0 * ff2_b
    end

    # Set and add follower moments initially resolved in basis b
    if has_mf_b
        # Set values
        @views mf1_b .= σ * mf_b[1,:,timeIndex]
        @views mf2_b .= σ * mf_b[2,:,timeIndex]
        # Add
        m1_ext .+= RR0 * mf1_b
        m2_ext .+= RR0 * mf2_b
    end

end

# Computes the nodal resultants from the aerodynamic loads and state arrays
function aerodynamic_loads_and_arrays!(f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,A,B,model,element,timeNow,Δt,isDynamicProblem)

    # Unpack states and rates
    @unpack p,V,Ω,χ = element.state.midpoint
    @unpack Vdot,Ωdot = element.state.rates
    @unpack δ = element.aero.surface.state

    # Unpack aerodynamic state cache
    @unpack stateCache = element.aero

    # Update aerodynamic nodal loads resutants and state arrays
    aero_core!(f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,A,B,model,element,stateCache,p,V,Ω,χ,δ,Vdot,Ωdot,timeNow,Δt,isDynamicProblem)

end

# Computes static contributions to the residual array for the current element
function element_static_residual!(F_u1,F_u2,F_p1,F_p2,F_F1,F_F2,F_M1,F_M2,Δℓ,R0,u,p,F,M,γ,κ,RR0,HTinv,f1,f2,m1,m2)

    # --- F_u --- #
    tmp = RR0 * F
    F_u1 .= -tmp - f1
    F_u2 .=  tmp - f2

    # --- F_p --- #
    tmp = RR0 * M
    tmp2 = Δℓ/2 * RR0 * cross(a1+γ, F)
    F_p1 .= -tmp - m1 - tmp2
    F_p2 .=  tmp - m2 - tmp2

    # --- F_F --- #
    tmp = Δℓ/2 * (RR0 * (a1+γ) - R0 * a1)
    F_F1 .=  u - tmp
    F_F2 .= -u - tmp

    # --- F_M --- #
    tmp = Δℓ/2 * HTinv * SVector(R0 * κ)
    F_M1 .=  p - tmp
    F_M2 .= -p - tmp

end

# Computes steady contributions to the residual array for the current element
function element_steady_residual!(F_u1,F_u2,F_p1,F_p2,F_V,F_Ω,F_χ,Δℓ,u,V,Ω,P,H,χ,v,ω,RR0,RR0T,A,B)

    # --- F_u --- #
    tmp = Δℓ/2 * cross(ω, RR0 * P)
    F_u1 .+= tmp
    F_u2 .+= tmp

    # --- F_p --- #
    tmp = Δℓ/2 * (cross(ω, RR0 * H) + RR0 * cross(V, P))
    F_p1 .+= tmp
    F_p2 .+= tmp

    # --- F_V --- #
    F_V .= RR0 * V - v - cross(ω, u)
    
    # --- F_Ω --- #
    F_Ω .= Ω - RR0T * ω

    # --- F_χ --- #
    F_χ .= - (A * χ + B)

end

# Computes transient contributions to the residual array for the current element
function element_transient_residual!(F_u1,F_u2,F_p1,F_p2,F_V,F_Ω,F_χ,Δℓ,R0T,P,H,udot,pdot,Pdot,Hdot,χdot,RR0,RdotR0,HT)

    # --- F_u --- #
    tmp = Δℓ/2 * (RdotR0 * P + RR0 * Pdot)
    F_u1 .+= tmp
    F_u2 .+= tmp

    # --- F_p --- #
    tmp = Δℓ/2 * (RdotR0 * H + RR0 * Hdot)
    F_p1 .+= tmp
    F_p2 .+= tmp

    # --- F_V --- #
    F_V .+= - udot
    
    # --- F_Ω --- #
    F_Ω .+= - R0T * SVector(HT * pdot)

    # --- F_χ --- #
    F_χ .+= χdot

end

# Inserts the element's contributions to the residual array
function insert_element_residual!(residual,F_u1,F_u2,F_p1,F_p2,F_F1,F_F2,F_M1,F_M2,F_V,F_Ω,F_χ,invForceScaling,eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,eqs_Fχ,isSpecialNode1,isSpecialNode2,eqsNode1Set,eqsNode2Set,hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat)

    # Convert to static vector for performance
    F_p1 = SVector(F_p1)

    # First node
    # --------------------------------------------------------------------------
    # The node's equations have not been set yet
    if !eqsNode1Set
        # Set equilibrium equations' residual
        residual[eqs_Fu1] .= F_u1 * invForceScaling
        residual[eqs_Fp1] .= F_p1 * invForceScaling
        # Set compatibility equations' residual
        residual[eqs_FF1] .= F_F1
        residual[eqs_FM1] .= F_M1
    # The node's equations have already been set
    else
        # Add to existing equilibrium equations' residual (except for hinged DOFs)
        residual[eqs_Fu1] .+= F_u1 * invForceScaling
        residual[eqs_Fp1] .+= notHingedNode1Mat * F_p1 * invForceScaling
        # This is a standard node (shares compatibility equations)
        if !isSpecialNode1
            # Add to existing compatibility equations' residual (except for hinged DOFs - discard equations) / set separate moment equilibrium equations' residual for hinged DOFs
            residual[eqs_FF1] .+= F_F1
            residual[eqs_FM1] .+= notHingedNode1Mat * F_M1 + 
                                  hingedNode1Mat * F_p1 * invForceScaling 
        # This is a special node (has separate compatibility equations)
        else
            # Set separate compatibility equations' residual (except for hinged DOFs - discard equations) / set separate moment equilibrium equations' residual for hinged DOFs
            residual[eqs_FF1] .= F_F1
            residual[eqs_FM1] .= notHingedNode1Mat * F_M1 + 
                                 hingedNode1Mat * F_p1 * invForceScaling
        end
    end

    # Second node
    # --------------------------------------------------------------------------
    # The node's equations have not been set yet
    if !eqsNode2Set
        # Set equilibrium equations' residual
        residual[eqs_Fu2] .= F_u2 * invForceScaling
        residual[eqs_Fp2] .= F_p2 * invForceScaling
        # Set compatibility equations' residual (except for hinged DOFs - discard equations)
        residual[eqs_FF2] .= F_F2
        residual[eqs_FM2] .= notHingedNode2Mat * F_M2 
    # The node's equations have already been set
    else
        # Add to existing equilibrium equations' residual
        residual[eqs_Fu2] .+= F_u2 * invForceScaling
        residual[eqs_Fp2] .+= F_p2 * invForceScaling
        # This is a standard node (shares compatibility equations)
        if !isSpecialNode2 
            # Add to existing compatibility equations' residual 
            residual[eqs_FF2] .+= F_F2
            residual[eqs_FM2] .+= F_M2
        # This is a special node (has separate compatibility equations)
        else
            # Set separate compatibility equations' residual 
            residual[eqs_FF2] .= F_F2
            residual[eqs_FM2] .= F_M2
        end
    end

    ## Midpoint
    # --------------------------------------------------------------------------
    # Generalized velocity-displacement
    residual[eqs_FV] .= F_V
    residual[eqs_FΩ] .= F_Ω

    # Aerodynamic
    residual[eqs_Fχ] .= F_χ

end

# Resets Jacobian entries that accumulate
function reset_jacobians!(jacobian,model,specialNodes,hingeAxisConstraints)

    # Reset spring loads' Jacobians to zero 
    reset_spring_loads_jacobians!(jacobian,specialNodes)

    # Reset trim flap deflections' columns to zero (all elements of a surface contribute to the same column)
    reset_trim_flap_deflections_jacobians!(jacobian,model)

end

# Resets cumulative Jacobians associated with spring loads
function reset_spring_loads_jacobians!(jacobian,specialNodes)

    @unpack nzval = jacobian

    # Loop special nodes
    for specialNode in specialNodes

        # Unpack special node data
        @unpack springs = specialNode
        @unpack doubleAttachmentSpringsData = specialNode.assemblyData
        @unpack Fu_uF,Fp_pM,Fu_uFother,Fp_pMother = specialNode.jacobianMaps

        # Initialize counter for double-attached springs
        c = 0
        
        # Loop special node's springs
        for spring in springs
            @unpack hasDoubleAttachment = spring.assemblyData
            if hasDoubleAttachment
                # Update counter
                c += 1
                # Reset Jacobians
                set_block!(nzval, Fu_uF, Z3)
                set_block!(nzval, Fp_pM, Z3)
                set_block!(nzval, Fu_uFother[c], Z3)
                set_block!(nzval, Fp_pMother[c], Z3)
            else
                set_block!(nzval, Fu_uF, Z3)
                set_block!(nzval, Fp_pM, Z3)
            end
        end

    end

end

# Resets cumulative Jacobians associated with the trim flap deflections
function reset_trim_flap_deflections_jacobians!(jacobian,model)

    @unpack colptr,nzval = jacobian
    @unpack rangeTrimδ = model.assemblyData

    # Loop trim flap deflections' columns
    for k in rangeTrimδ
        k < length(colptr) || break
        @views nzval[colptr[k]:colptr[k+1]-1] .= 0
    end

end

# Computes the contributions from the current element to the Jacobian matrix
function element_jacobian!(jacobian,model,element,timeNow,Δt,isDynamicProblem,skipAeroJacobianUpdate)

    @unpack forceScaling,invForceScaling = model.assemblyData
    @unpack jacobianMaps = element
    @unpack Δℓ,R0,R0T,hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat = element.geometry
    @unpack C11,C12,C21,C22,I11,I12,I21,I22 = element.sectional
    @unpack F,M,V,Ω,γ,κ,P,H,χ = element.state.midpoint
    @unpack pdot,Vdot,Ωdot,Pdot,Hdot,χdot = element.state.rates
    @unpack ω = element.stateCache.kinematicsBasisb
    @unpack RR0,RdotR0,HT,HTinv,R_p1,R_p2,R_p3,RT_p1,RT_p2,RT_p3,HT_p1,HT_p2,HT_p3,HTinv_p1,HTinv_p2,HTinv_p3,Rdot_p1,Rdot_p2,Rdot_p3 = element.stateCache.rotationMatrices
    @unpack f1_p,f2_p,m1_p,m2_p = element.stateCache.loadsJacobians
    @unpack f1χ_V,f2χ_V,f1χ_Ω,f2χ_Ω,m1χ_V,m2χ_V,m1χ_Ω,m2χ_Ω,f1χ_χ,f2χ_χ,m1χ_χ,m2χ_χ,f1χ_δ,f2χ_δ,m1χ_δ,m2χ_δ,f1χ_Vdot,f2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Vdot,m2χ_Vdot,m1χ_Ωdot,m2χ_Ωdot = element.stateCache.aeroLoadsJacobians
    @unpack F_u1_p,F_u2_p,F_u1_F,F_u2_F,F_u1_V,F_u2_V,F_u1_Ω,F_u2_Ω,F_u1_χ,F_u2_χ,F_u1_δ,F_u2_δ,F_p1_p,F_p2_p,F_p1_F,F_p2_F,F_p1_M,F_p2_M,F_p1_V,F_p2_V,F_p1_Ω,F_p2_Ω,F_p1_χ,F_p2_χ,F_p1_δ,F_p2_δ,F_F1_u,F_F2_u,F_F1_p,F_F2_p,F_F1_F,F_F2_F,F_F1_M,F_F2_M,F_M1_p,F_M2_p,F_M1_F,F_M2_F,F_M1_M,F_M2_M,F_V_u,F_V_p,F_V_V,F_Ω_p,F_Ω_Ω,F_χ_p,F_χ_V,F_χ_Ω,F_χ_χ,F_χ_Vdot,F_χ_Ωdot = element.stateCache.jacobians

    # Derivatives of the rotation variables w.r.t. the rotation parameters
    element_rotation_derivatives!(element,Δt,isDynamicProblem)

    # Convert to static matrices for performance
    RR0 = SMatrix(RR0)

    # Aerodynamic derivatives (contributions to Jacobian matrix)
    if !skipAeroJacobianUpdate 
        aero_derivatives_states!(model,element,timeNow,Δt,isDynamicProblem)
        isDynamicProblem && aero_derivatives_states_rates!(model,element,timeNow,Δt,isDynamicProblem)
    end

    # Follower distributed loads' derivatives (including aerodynamic loads) w.r.t. the rotation parameters (contributions to Jacobian matrix)
    distributed_loads_derivatives_rotation_parameters!(element)

    # Static terms
    element_static_jacobians!(F_u1_p,F_u2_p,F_u1_F,F_u2_F,F_p1_p,F_p2_p,F_p1_F,F_p2_F,F_p1_M,F_p2_M,F_F1_u,F_F2_u,F_F1_p,F_F2_p,F_F1_F,F_F2_F,F_F1_M,F_F2_M,F_M1_p,F_M2_p,F_M1_F,F_M2_F,F_M1_M,F_M2_M,Δℓ,R0,C11,C12,C21,C22,F,M,γ,κ,RR0,HTinv,R_p1,R_p2,R_p3,HTinv_p1,HTinv_p2,HTinv_p3,f1_p,f2_p,m1_p,m2_p)

    # Steady-state terms
    element_steady_jacobians!(F_u1_p,F_u2_p,F_u1_V,F_u2_V,F_u1_Ω,F_u2_Ω,F_u1_χ,F_u2_χ,F_u1_δ,F_u2_δ,F_p1_p,F_p2_p,F_p1_V,F_p2_V,F_p1_Ω,F_p2_Ω,F_p1_χ,F_p2_χ,F_p1_δ,F_p2_δ,F_V_u,F_V_p,F_V_V,F_Ω_p,F_Ω_Ω,Δℓ,I11,I12,I21,I22,R0,R0T,ω,V,P,H,RR0,R_p1,R_p2,R_p3,RT_p1,RT_p2,RT_p3,f1χ_V,f2χ_V,f1χ_Ω,f2χ_Ω,m1χ_V,m2χ_V,m1χ_Ω,m2χ_Ω,f1χ_χ,f2χ_χ,m1χ_χ,m2χ_χ,f1χ_δ,f2χ_δ,m1χ_δ,m2χ_δ)

    # Transient dynamic terms
    isDynamicProblem && element_transient_jacobians!(F_u1_p,F_u2_p,F_u1_V,F_u2_V,F_u1_Ω,F_u2_Ω,F_p1_p,F_p2_p,F_p1_V,F_p2_V,F_p1_Ω,F_p2_Ω,F_V_u,F_Ω_p,F_χ_V,F_χ_Ω,F_χ_χ,F_χ_Vdot,F_χ_Ωdot,f1χ_Vdot,f2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Vdot,m2χ_Vdot,m1χ_Ωdot,m2χ_Ωdot,Δt,Δℓ,I11,I12,I21,I22,R0,R0T,P,H,pdot,Pdot,Hdot,RR0,RdotR0,HT,R_p1,R_p2,R_p3,HT_p1,HT_p2,HT_p3,Rdot_p1,Rdot_p2,Rdot_p3)

    # Insert element resultants into the Jacobian matrix
    insert_element_jacobians!(jacobian,F_u1_p,F_u2_p,F_u1_F,F_u2_F,F_u1_V,F_u2_V,F_u1_Ω,F_u2_Ω,F_u1_χ,F_u2_χ,F_u1_δ,F_u2_δ,F_p1_p,F_p2_p,F_p1_F,F_p2_F,F_p1_M,F_p2_M,F_p1_V,F_p2_V,F_p1_Ω,F_p2_Ω,F_p1_χ,F_p2_χ,F_p1_δ,F_p2_δ,F_F1_u,F_F2_u,F_F1_p,F_F2_p,F_F1_F,F_F2_F,F_F1_M,F_F2_M,F_M1_p,F_M2_p,F_M1_F,F_M2_F,F_M1_M,F_M2_M,F_V_u,F_V_p,F_V_V,F_Ω_p,F_Ω_Ω,F_χ_p,F_χ_V,F_χ_Ω,F_χ_χ,forceScaling,invForceScaling,hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat,jacobianMaps)

    # @logmsg LogLevel(show_debug() ? 0 : -1) "ELEMENT $(element.assemblyData.globalID) JACOBIANS"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_u1_p: $F_u1_p"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_u1_F: $F_u1_F"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_u1_V: $F_u1_V"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_u1_Ω: $F_u1_Ω"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p1_p: $F_p1_p"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p1_F: $F_p1_F"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p1_M: $F_p1_M"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p1_V: $F_p1_V"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p1_Ω: $F_p1_Ω"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_F1_u: $F_F1_u"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_F1_p: $F_F1_p"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_F1_F: $F_F1_F"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_F1_M: $F_F1_M"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_M1_p: $F_M1_p"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_M1_F: $F_M1_F"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_M1_M: $F_M1_M"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_V_u: $F_V_u"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_V_p: $F_V_p"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_V_V: $F_V_V"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_Ω_p: $F_Ω_p"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_Ω_Ω: $F_Ω_Ω"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_χ_p: $F_χ_p"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_χ_V: $F_χ_V"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_χ_Ω: $F_χ_Ω"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_χ_χ: $F_χ_χ"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_u1_δ: $F_u1_δ"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_u2_δ: $F_u2_δ"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p1_δ: $F_p1_δ"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p2_δ: $F_p2_δ"
    # @logmsg LogLevel(show_debug() ? 0 : -1) (norm(jacobian))
    # if size(jacobian, 1) == size(jacobian, 2)
    #     @logmsg LogLevel(show_debug() ? 0 : -1) (det(jacobian))
    # end

    return nothing
    
end

# Computes the derivatives of the rotation variables w.r.t. the extended rotation parameters (required by the Jacobian and inertia matrices only)
function element_rotation_derivatives!(element,Δt,isDynamicProblem)

    @unpack p = element.state.midpoint
    @unpack pdot = element.state.rates
    @unpack R,HT,R_p1,R_p2,R_p3,RT_p1,RT_p2,RT_p3,HT_p1,HT_p2,HT_p3,HTinv_p1,HTinv_p2,HTinv_p3,Rdot_p1,Rdot_p2,Rdot_p3 = element.stateCache.rotationMatrices

    # Derivatives of the rotation tensor and the tangent operator tensors w.r.t. extended rotation parameters
    ∂R_∂p_and_∂HT_∂p_and_∂HTinv_∂p!(R_p1,R_p2,R_p3,HT_p1,HT_p2,HT_p3,HTinv_p1,HTinv_p2,HTinv_p3,p)

    # Transposes
    RT_p1 .= R_p1'
    RT_p2 .= R_p2'
    RT_p3 .= R_p3'

    # Derivatives of the time derivative of rotation tensor for dynamic analyses
    isDynamicProblem && ∂Rdot_∂p!(Rdot_p1,Rdot_p2,Rdot_p3,R,HT,R_p1,R_p2,R_p3,HT_p1,HT_p2,HT_p3,pdot,Δt)

end

# Computes the derivatives of the distributed loads w.r.t the rotation parameters
function distributed_loads_derivatives_rotation_parameters!(element)

    @unpack f1_ext_p,f2_ext_p,m1_ext_p,m2_ext_p,mg_p = element.stateCache.structuralLoadsJacobians
    @unpack f1χ_p,f2χ_p,m1χ_p,m2χ_p = element.stateCache.aeroLoadsJacobians
    @unpack f1_p,f2_p,m1_p,m2_p = element.stateCache.loadsJacobians

    # Derivatives of gravitational loads w.r.t. extended rotation parameters
    gravitational_loads_derivatives_rotation_parameters!(mg_p,element)

    # Derivatives of externally applied follower distributed loads w.r.t. extended rotation parameters
    distributed_external_loads_derivatives_rotation_parameters!(f1_ext_p,f2_ext_p,m1_ext_p,m2_ext_p,element)

    # Derivatives of aerodynamic loads w.r.t. extended rotation parameters (f1χ_p,f2χ_p,m1χ_p,m2χ_p) already computed by aero_derivatives_states!    

    # Total nodal resultants' derivatives w.r.t. extended rotation parameters
    f1_p .= f1_ext_p + f1χ_p
    f2_p .= f2_ext_p + f2χ_p
    m1_p .= m1_ext_p + m1χ_p + mg_p
    m2_p .= m2_ext_p + m2χ_p + mg_p

end


# Computes the contributions of the gravitational loads to the Jacobian matrix
function gravitational_loads_derivatives_rotation_parameters!(mg_p,element)

    @unpack R0T = element.geometry
    @unpack ηtilde = element.sectional
    @unpack R,RR0,R_p1,R_p2,R_p3,RT_p1,RT_p2,RT_p3 = element.stateCache.rotationMatrices
    @unpack f_g,m_g = element.stateCache.structuralLoadsResultants

    # Derivatives of gravitational loads w.r.t. extended rotation parameters
    if any(!iszero(f_g))
        mg_p .= mul3(R_p1,R_p2,R_p3,SVector(R'*m_g)) + 
            SMatrix(RR0 * ηtilde) * SMatrix(R0T * mul3(RT_p1,RT_p2,RT_p3,f_g))
    end

end


# Computes the contributions of the externally applied distributed loads to the Jacobian matrix
function distributed_external_loads_derivatives_rotation_parameters!(f1_ext_p,f2_ext_p,m1_ext_p,m2_ext_p,element)

    @unpack R0 = element.geometry
    @unpack R_p1,R_p2,R_p3 = element.stateCache.rotationMatrices
    @unpack has_ff_A,has_mf_A,has_ff_b,has_mf_b = element.distLoads
    @unpack ff1_A,ff2_A,mf1_A,mf2_A,ff1_b,ff2_b,mf1_b,mf2_b = element.stateCache.structuralLoadsResultants
    
    # Reset to zero arrays
    zero_arrays!(f1_ext_p,f2_ext_p,m1_ext_p,m2_ext_p)

    # Contributions from follower forces initially resolved in basis A
    if has_ff_A
        f1_ext_p .+= mul3(R_p1,R_p2,R_p3,ff1_A)
        f2_ext_p .+= mul3(R_p1,R_p2,R_p3,ff2_A)
    end

    # Contributions from follower moments initially resolved in basis A
    if has_mf_A
        m1_ext_p .+= mul3(R_p1,R_p2,R_p3,mf1_A)
        m2_ext_p .+= mul3(R_p1,R_p2,R_p3,mf2_A)
    end

    # Contributions from follower forces initially resolved in basis b
    if has_ff_b
        f1_ext_p .+= mul3(R_p1,R_p2,R_p3,R0*ff1_b)
        f2_ext_p .+= mul3(R_p1,R_p2,R_p3,R0*ff2_b)
    end

    # Contributions from follower moements initially resolved in basis b
    if has_mf_b
        m1_ext_p .+= mul3(R_p1,R_p2,R_p3,R0*mf1_b)
        m2_ext_p .+= mul3(R_p1,R_p2,R_p3,R0*mf2_b)
    end

end

# Computes the derivatives of the aerodynamic loads and aerodynamic states w.r.t elemental states (excludes rotation parameters, already computed)
function aero_derivatives_states!(model,element,timeNow,Δt,isDynamicProblem)

    @unpack aero = element

    # Skip if there is no aerodynamic surface, or angle of attack is undefined
    if isnothing(aero.surface) || isnan(aero.stateCache.attachedVars.steadyKin.α)
        @unpack f1χ_p,f2χ_p,f1χ_V,f2χ_V,f1χ_Ω,f2χ_Ω,m1χ_p,m2χ_p,m1χ_V,m2χ_V,m1χ_Ω,m2χ_Ω,f1χ_χ,f2χ_χ,m1χ_χ,m2χ_χ,f1χ_δ,f2χ_δ,m1χ_δ,m2χ_δ = element.stateCache.aeroLoadsJacobians
        @unpack F_χ_p,F_χ_V,F_χ_Ω = element.stateCache.jacobians
        zero_arrays!(f1χ_p,f2χ_p,f1χ_V,f2χ_V,f1χ_Ω,f2χ_Ω,m1χ_p,m2χ_p,m1χ_V,m2χ_V,m1χ_Ω,m2χ_Ω,f1χ_χ,f2χ_χ,m1χ_χ,m2χ_χ,f1χ_δ,f2χ_δ,m1χ_δ,m2χ_δ,F_χ_p,F_χ_V,F_χ_Ω)
        return
    end

    # Unpack element data
    @unpack p,V,Ω,χ = element.state.midpoint
    @unpack δ = element.aero.surface.state
    @unpack derivationMethod = aero.surface.solvers
    @unpack DOF_δ = aero.surface.assemblyData

    # Input states for aero loads and matrices wrapper function
    states = aero.stateCache.statesBuffer
    copyto!(view(states,1:3), p)
    copyto!(view(states,4:6), V)
    copyto!(view(states,7:9), Ω)
    copyto!(view(states,10:9+length(χ)), χ)
    if !isempty(DOF_δ)
        @inbounds states[end] = δ
    end

    # Jacobians of aerodynamic loads and state matrices w.r.t. states
    aero_jacobians_states!(states,model,element,timeNow,Δt,isDynamicProblem,derivationMethod)

end

# Computes the aerodynamic Jacobians of states with Automatic Differentiation
function aero_jacobians_states!(states,model,element,timeNow,Δt,isDynamicProblem,derivationMethod::AD)

    # Derivatives of aerodynamic loads and state matrices w.r.t. states
    J = ForwardDiff.jacobian(x -> wrapper_aero_loads_and_matrices_from_states!(x,model,element,timeNow,Δt,isDynamicProblem), SVector(states))

    # Set derivatives
    set_aero_jacobians_of_states!(J,element)

end

# Computes the aerodynamic Jacobians of states with Finite Differences
function aero_jacobians_states!(states,model,element,timeNow,Δt,isDynamicProblem,derivationMethod::FD)

    # Derivatives of aerodynamic loads and state matrices w.r.t. states
    J = first(FiniteDifferences.jacobian(derivationMethod.method, x -> wrapper_aero_loads_and_matrices_from_states!(x,model,element,timeNow,Δt,isDynamicProblem), states))

    # Set derivatives
    set_aero_jacobians_of_states!(J,element)

end

# Wrapper function for aerodynamic loads and state matrices using elemental states as inputs
function wrapper_aero_loads_and_matrices_from_states!(states,model,element,timeNow,Δt,isDynamicProblem)

    # Unpack states and rates
    p = @view states[1:3]
    V = @view states[4:6]
    Ω = @view states[7:9]
    χ = @view states[10:9+length(element.state.midpoint.χ)]
    δ = !isempty(element.aero.surface.assemblyData.DOF_δ) ? 
        states[end] :
        element.aero.surface.state.δ
    @unpack Vdot,Ωdot = element.state.rates

    # Set state cache for derivatives calculations
    stateCacheDiff = get_diff_state_cache!(element.aero,eltype(states))

    # Unpack resultants and arrays
    @unpack f1χ,f2χ,m1χ,m2χ = stateCacheDiff.loadsResultants
    @unpack Ac,Ai,An,Af,Ag,A,B = stateCacheDiff.stateMatrices

    # Compute nodal resultants and state arrays
    f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,B = aero_core!(f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,A,B,model,element,stateCacheDiff,p,V,Ω,χ,δ,Vdot,Ωdot,timeNow,Δt,isDynamicProblem)

    # Set outputs array
    outBuffer = stateCacheDiff.outBuffer
    write_aero_jacobians_outputs!(outBuffer, f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,B)

    return outBuffer

end

# Writes the aerodynamic loads and state matrices into the output buffer
@inline function write_aero_jacobians_outputs!(out, f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,B)
    @inbounds begin
        copyto!(view(out, 1:3),   f1χ)
        copyto!(view(out, 4:6),   f2χ)
        copyto!(view(out, 7:9),   m1χ)
        copyto!(view(out, 10:12), m2χ)

        i = 13

        n = length(Ac); copyto!(view(out, i:i+n-1), Ac); i += n
        n = length(Ai); copyto!(view(out, i:i+n-1), Ai); i += n
        n = length(An); copyto!(view(out, i:i+n-1), An); i += n
        n = length(Af); copyto!(view(out, i:i+n-1), Af); i += n
        n = length(Ag); copyto!(view(out, i:i+n-1), Ag); i += n
        n = length(B);  copyto!(view(out, i:i+n-1), B)
    end
    return out
end

# Computes the aerodynamic derivatives of the element's states
function set_aero_jacobians_of_states!(J,element)

    # Unpack element data
    @unpack χ = element.state.midpoint
    @unpack ppSolver = element.aero.surface.solvers
    @unpack DOF_δ = element.aero.surface.assemblyData
    @unpack f1χ_p,f2χ_p,f1χ_V,f2χ_V,f1χ_Ω,f2χ_Ω,m1χ_p,m2χ_p,m1χ_V,m2χ_V,m1χ_Ω,m2χ_Ω,f1χ_χ,f2χ_χ,m1χ_χ,m2χ_χ,f1χ_δ,f2χ_δ,m1χ_δ,m2χ_δ = element.stateCache.aeroLoadsJacobians
    @unpack F_χ_p,F_χ_V,F_χ_Ω,F_χ_χ = element.stateCache.jacobians
    @unpack A,B,Ac,Ai,An,Af,Ag = element.aero.stateCache.stateMatrices
    @unpack Ac_p,Ai_p,An_p,Af_p,Ag_p,Ac_V,Ai_V,An_V,Af_V,Ag_V,Ac_Ω,Ai_Ω,An_Ω,Af_Ω,Ag_Ω,Ac_χ,Ai_χ,An_χ,Af_χ,Ag_χ,A_p,A_V,A_Ω,A_χ,B_p,B_V,B_Ω,B_χ = element.aero.matricesDerivatives

    # Derivatives of aerodynamic nodal loads resultants w.r.t. p
    f1χ_p .= @view J[1:3, 1:3]
    f2χ_p .= @view J[4:6, 1:3]
    m1χ_p .= @view J[7:9, 1:3]
    m2χ_p .= @view J[10:12, 1:3]

    # Derivatives of aerodynamic nodal loads resultants w.r.t. V
    f1χ_V .= @view J[1:3, 4:6]
    f2χ_V .= @view J[4:6, 4:6]
    m1χ_V .= @view J[7:9, 4:6]
    m2χ_V .= @view J[10:12, 4:6]

    # Derivatives of aerodynamic nodal loads resultants w.r.t. Ω
    f1χ_Ω .= @view J[1:3, 7:9]
    f2χ_Ω .= @view J[4:6, 7:9]
    m1χ_Ω .= @view J[7:9, 7:9]
    m2χ_Ω .= @view J[10:12, 7:9]

    # Derivatives of aerodynamic nodal loads resultants w.r.t. χ
    if !isempty(χ)
        f1χ_χ .= @view J[1:3, 10:9+length(χ)]
        f2χ_χ .= @view J[4:6, 10:9+length(χ)]
        m1χ_χ .= @view J[7:9, 10:9+length(χ)]
        m2χ_χ .= @view J[10:12, 10:9+length(χ)]
    end

    # Derivatives of aerodynamic nodal loads resultants w.r.t. δ
    if !isempty(DOF_δ)
        f1χ_δ .= @view J[1:3, end]
        f2χ_δ .= @view J[4:6, end]
        m1χ_δ .= @view J[7:9, end]
        m2χ_δ .= @view J[10:12, end]
    end

    # Derivatives of aerodynamic equations w.r.t. V and Ω
    if !isempty(χ)

        # Lengths of arrays
        Nc = length(Ac)
        Ni = length(Ai)
        Nn = length(An)
        Nf = length(Af)
        Ng = length(Ag)

        # Derivatives of state arrays
        Ac_p .= @view J[13 : 12+Nc, 1:3]
        Ac_V .= @view J[13 : 12+Nc, 4:6]
        Ac_Ω .= @view J[13 : 12+Nc, 7:9]
        Ac_χ .= @view J[13 : 12+Nc, 10:9+length(χ)]
        Ai_p .= @view J[13+Nc : 12+Nc+Ni, 1:3]
        Ai_V .= @view J[13+Nc : 12+Nc+Ni, 4:6]
        Ai_Ω .= @view J[13+Nc : 12+Nc+Ni, 7:9]
        Ai_χ .= @view J[13+Nc : 12+Nc+Ni, 10:9+length(χ)]
        An_p .= @view J[13+Nc+Ni : 12+Nc+Ni+Nn, 1:3]
        An_V .= @view J[13+Nc+Ni : 12+Nc+Ni+Nn, 4:6]
        An_Ω .= @view J[13+Nc+Ni : 12+Nc+Ni+Nn, 7:9]
        An_χ .= @view J[13+Nc+Ni : 12+Nc+Ni+Nn, 10:9+length(χ)]
        Af_p .= @view J[13+Nc+Ni+Nn : 12+Nc+Ni+Nn+Nf, 1:3]
        Af_V .= @view J[13+Nc+Ni+Nn : 12+Nc+Ni+Nn+Nf, 4:6]
        Af_Ω .= @view J[13+Nc+Ni+Nn : 12+Nc+Ni+Nn+Nf, 7:9]
        Af_χ .= @view J[13+Nc+Ni+Nn : 12+Nc+Ni+Nn+Nf, 10:9+length(χ)]
        Ag_p .= @view J[13+Nc+Ni+Nn+Nf : 12+Nc+Ni+Nn+Nf+Ng, 1:3]
        Ag_V .= @view J[13+Nc+Ni+Nn+Nf : 12+Nc+Ni+Nn+Nf+Ng, 4:6]
        Ag_Ω .= @view J[13+Nc+Ni+Nn+Nf : 12+Nc+Ni+Nn+Nf+Ng, 7:9]        
        Ag_χ .= @view J[13+Nc+Ni+Nn+Nf : 12+Nc+Ni+Nn+Nf+Ng, 10:9+length(χ)]

        # Set A_p, A_V, A_Ω and A_χ
        @inbounds for k = 1:3
            @views begin
                i = 1
                if ppSolver isa Inflow
                    n = isqrt(Nc)
                    A_p[i:i+n-1, i:i+n-1, k] .= reshape(Ac_p[:,k], (n,n))
                    A_V[i:i+n-1, i:i+n-1, k] .= reshape(Ac_V[:,k], (n,n))
                    A_Ω[i:i+n-1, i:i+n-1, k] .= reshape(Ac_Ω[:,k], (n,n))
                    i += n
                else
                    i = assemble_A_diagonal_block!((A_p,A_V,A_Ω), (Ac_p[:,k],Ac_V[:,k],Ac_Ω[:,k]), i, k)
                end
                i = assemble_A_diagonal_block!((A_p,A_V,A_Ω), (Ai_p[:,k],Ai_V[:,k],Ai_Ω[:,k]), i, k)
                i = assemble_A_diagonal_block!((A_p,A_V,A_Ω), (An_p[:,k],An_V[:,k],An_Ω[:,k]), i, k)
                i = assemble_A_diagonal_block!((A_p,A_V,A_Ω), (Af_p[:,k],Af_V[:,k],Af_Ω[:,k]), i, k)
                i = assemble_A_diagonal_block!((A_p,A_V,A_Ω), (Ag_p[:,k],Ag_V[:,k],Ag_Ω[:,k]), i, k)
            end
        end
        @inbounds for k in eachindex(χ)
            @views begin
                i = 1
                if ppSolver isa Inflow
                    n = isqrt(Nc)
                    A_χ[i:i+n-1, i:i+n-1, k] .= reshape(Ac_χ[:,k], (n,n))
                    i += n
                else
                    i = assemble_A_diagonal_block!((A_χ,), (Ac_χ[:,k],), i, k)
                end
                i = assemble_A_diagonal_block!((A_χ,), (Ai_χ[:,k],), i, k)
                i = assemble_A_diagonal_block!((A_χ,), (An_χ[:,k],), i, k)
                i = assemble_A_diagonal_block!((A_χ,), (Af_χ[:,k],), i, k)
                i = assemble_A_diagonal_block!((A_χ,), (Ag_χ[:,k],), i, k)
            end
        end

        # Set B_p, B_V, B_Ω and B_χ
        B_p .= @view J[end-length(χ)+1:end, 1:3]
        B_V .= @view J[end-length(χ)+1:end, 4:6]
        B_Ω .= @view J[end-length(χ)+1:end, 7:9]
        B_χ .= @view J[end-length(χ)+1:end, 10:9+length(χ)]

        # Set Jacobians
        @inbounds for i=1:3
            F_χ_p[:,i] .= - (A_p[:,:,i] * χ + B_p[:,i])
            F_χ_V[:,i] .= - (A_V[:,:,i] * χ + B_V[:,i])
            F_χ_Ω[:,i] .= - (A_Ω[:,:,i] * χ + B_Ω[:,i])
        end
        @inbounds for i in eachindex(χ)
            F_χ_χ[:,i] .= -(A[:,i] + A_χ[:,:,i] * χ + B_χ[:,i])
        end
    end
    
end

# Assembles a diagonal block into an arbitrary number of derivative arrays A_x[:,:,k]
@inline function assemble_A_diagonal_block!(As::Tuple, vs::Tuple, i, k)
    n = length(first(vs))
    @inbounds for j = 1:n
        for m in eachindex(As)
            As[m][i+j-1, i+j-1, k] = vs[m][j]
        end
    end
    return i + n
end

# Computes static contributions to the Jacobian matrix from the current element
function element_static_jacobians!(F_u1_p,F_u2_p,F_u1_F,F_u2_F,F_p1_p,F_p2_p,F_p1_F,F_p2_F,F_p1_M,F_p2_M,F_F1_u,F_F2_u,F_F1_p,F_F2_p,F_F1_F,F_F2_F,F_F1_M,F_F2_M,F_M1_p,F_M2_p,F_M1_F,F_M2_F,F_M1_M,F_M2_M,Δℓ,R0,C11,C12,C21,C22,F,M,γ,κ,RR0,HTinv,R_p1,R_p2,R_p3,HTinv_p1,HTinv_p2,HTinv_p3,f1_p,f2_p,m1_p,m2_p)

    # Repeatedly used variables
    F̃ = tilde(F)
        
    # --- F_u --- #
    # F_u_p
    tmp = mul3(R_p1,R_p2,R_p3,SVector(R0*F))
    F_u1_p .= -tmp - f1_p
    F_u2_p .=  tmp - f2_p
    # F_u_F
    F_u1_F .= -RR0
    F_u2_F .=  RR0

    # --- F_p --- #
    # F_p_p
    tmp1 = mul3(R_p1,R_p2,R_p3,SVector(R0*M))
    tmp2 = Δℓ/2 * mul3(R_p1,R_p2,R_p3,R0*cross(a1+γ,F))
    F_p1_p .= -tmp1 - m1_p - tmp2
    F_p2_p .=  tmp1 - m2_p - tmp2
    # F_p_F
    F_p1_F .= -Δℓ/2 * RR0 * (tilde(a1+γ) - F̃ * C11)
    F_p2_F .= F_p1_F
    # F_p_M
    tmp = Δℓ/2 * RR0 * F̃ * C12
    F_p1_M .= tmp - RR0
    F_p2_M .= tmp + RR0

    # --- F_F --- #
    # F_F_u (already set as constants)
    # F_F1_u .=  I3
    # F_F2_u .= -I3
    # F_F_p
    F_F1_p .= -Δℓ/2 * mul3(R_p1,R_p2,R_p3,R0*(a1+γ))
    F_F2_p .= F_F1_p
    # F_F_F
    F_F1_F .= -Δℓ/2 * RR0 * C11
    F_F2_F .= F_F1_F
    # F_F_M
    F_F1_M .= -Δℓ/2 * RR0 * C12
    F_F2_M .= F_F1_M

    # --- F_M --- #
    # F_M_p
    tmp = Δℓ/2 * mul3(HTinv_p1,HTinv_p2,HTinv_p3,SVector(R0*κ))
    F_M1_p .=  I3 - tmp
    F_M2_p .= -I3 - tmp
    # F_M_F
    tmp = -Δℓ/2 * SMatrix(HTinv * R0)
    F_M1_F .= tmp * C21
    F_M2_F .= F_M1_F
    # F_M_M
    F_M1_M .= tmp * C22
    F_M2_M .= F_M1_M

end

# Computes steady-state contributions to the Jacobian matrix from the current element
function element_steady_jacobians!(F_u1_p,F_u2_p,F_u1_V,F_u2_V,F_u1_Ω,F_u2_Ω,F_u1_χ,F_u2_χ,F_u1_δ,F_u2_δ,F_p1_p,F_p2_p,F_p1_V,F_p2_V,F_p1_Ω,F_p2_Ω,F_p1_χ,F_p2_χ,F_p1_δ,F_p2_δ,F_V_u,F_V_p,F_V_V,F_Ω_p,F_Ω_Ω,Δℓ,I11,I12,I21,I22,R0,R0T,ω,V,P,H,RR0,R_p1,R_p2,R_p3,RT_p1,RT_p2,RT_p3,f1χ_V,f2χ_V,f1χ_Ω,f2χ_Ω,m1χ_V,m2χ_V,m1χ_Ω,m2χ_Ω,f1χ_χ,f2χ_χ,m1χ_χ,m2χ_χ,f1χ_δ,f2χ_δ,m1χ_δ,m2χ_δ)

    # Repeatedly used variables
    ω̃ = tilde(ω)
    ω̃RR0 = ω̃ * RR0
    Ṽ = tilde(V)

    # --- F_u --- #
    # F_u_p
    tmp = Δℓ/2 * ω̃ * mul3(R_p1,R_p2,R_p3,SVector(R0*P))
    F_u1_p .+= tmp
    F_u2_p .+= tmp
    # F_u_V
    tmp = Δℓ/2 * ω̃RR0 * I11
    F_u1_V .= tmp - f1χ_V
    F_u2_V .= tmp - f2χ_V
    # F_u_Ω
    tmp = Δℓ/2 * ω̃RR0 * I12
    F_u1_Ω .= tmp - f1χ_Ω
    F_u2_Ω .= tmp - f2χ_Ω
    # F_u_χ
    F_u1_χ .= - f1χ_χ
    F_u2_χ .= - f2χ_χ
    # F_u_δ
    F_u1_δ .= - f1χ_δ
    F_u2_δ .= - f2χ_δ
    
    # --- F_p --- #
    # F_p_p
    tmp = Δℓ/2 * (ω̃ * mul3(R_p1,R_p2,R_p3,SVector(R0*H)) + 
                  mul3(R_p1,R_p2,R_p3,SVector(R0*cross(V,P))))
    F_p1_p .+= tmp 
    F_p2_p .+= tmp
    # F_p_V
    tmp = Δℓ/2 * (ω̃RR0 * I21 + RR0 * ( Ṽ * I11 - tilde(P) ))
    F_p1_V .= tmp - m1χ_V
    F_p2_V .= tmp - m2χ_V
    # F_p_Ω
    tmp = Δℓ/2 * (ω̃RR0 * I22 + RR0 * Ṽ * I12)
    F_p1_Ω .= tmp - m1χ_Ω
    F_p2_Ω .= tmp - m2χ_Ω
    # F_p_χ
    F_p1_χ .= - m1χ_χ
    F_p2_χ .= - m2χ_χ
    # F_p_δ
    F_p1_δ .= - m1χ_δ
    F_p2_δ .= - m2χ_δ
    
    # --- F_V --- #
    # F_V_u
    F_V_u .= - ω̃
    # F_V_p
    F_V_p .= mul3(R_p1,R_p2,R_p3,SVector(R0*V))
    # F_V_V
    F_V_V .= RR0
    
    # --- F_Ω --- #
    # F_Ω_p
    F_Ω_p .= - R0T * mul3(RT_p1,RT_p2,RT_p3,ω)
    # F_Ω_Ω (already set as constant)
    # F_Ω_Ω .= I3

end

# Computes transient dynamic contributions to the Jacobian matrix from the current element
function element_transient_jacobians!(F_u1_p,F_u2_p,F_u1_V,F_u2_V,F_u1_Ω,F_u2_Ω,F_p1_p,F_p2_p,F_p1_V,F_p2_V,F_p1_Ω,F_p2_Ω,F_V_u,F_Ω_p,F_χ_V,F_χ_Ω,F_χ_χ,F_χ_Vdot,F_χ_Ωdot,f1χ_Vdot,f2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Vdot,m2χ_Vdot,m1χ_Ωdot,m2χ_Ωdot,Δt,Δℓ,I11,I12,I21,I22,R0,R0T,P,H,pdot,Pdot,Hdot,RR0,RdotR0,HT,R_p1,R_p2,R_p3,HT_p1,HT_p2,HT_p3,Rdot_p1,Rdot_p2,Rdot_p3)

    # Repeatedly used variables
    RdotR0_plus_2oΔtRR0 = SMatrix(RdotR0 + 2/Δt*RR0)
    
    # --- F_u --- #
    # F_u_p
    tmp = Δℓ/2 * (mul3(Rdot_p1,Rdot_p2,Rdot_p3,SVector(R0*P)) + 
                  mul3(R_p1,R_p2,R_p3,SVector(R0*Pdot)))
    F_u1_p .+= tmp
    F_u2_p .+= tmp
    # F_u_V
    tmp = Δℓ/2 * RdotR0_plus_2oΔtRR0 * I11
    F_u1_V .+= tmp - f1χ_Vdot * 2/Δt
    F_u2_V .+= tmp - f2χ_Vdot * 2/Δt
    # F_u_Ω
    tmp = Δℓ/2 * RdotR0_plus_2oΔtRR0 * I12
    F_u1_Ω .+= tmp - f1χ_Ωdot * 2/Δt
    F_u2_Ω .+= tmp - f2χ_Ωdot * 2/Δt

    # --- F_p --- #
    # F_p_p
    tmp = Δℓ/2 * (mul3(Rdot_p1,Rdot_p2,Rdot_p3,SVector(R0*H)) + 
                  mul3(R_p1,R_p2,R_p3,SVector(R0*Hdot)))
    F_p1_p .+= tmp
    F_p2_p .+= tmp
    # F_p_V
    tmp = Δℓ/2 * RdotR0_plus_2oΔtRR0 * I21
    F_p1_V .+= tmp - m1χ_Vdot * 2/Δt
    F_p2_V .+= tmp - m2χ_Vdot * 2/Δt
    # F_p_Ω
    tmp = Δℓ/2 * RdotR0_plus_2oΔtRR0 * I22
    F_p1_Ω .+= tmp - m1χ_Ωdot * 2/Δt
    F_p2_Ω .+= tmp - m2χ_Ωdot * 2/Δt
    
    # --- F_V --- #
    # F_V_u
    F_V_u .+= - 2/Δt * I3

    # --- F_Ω --- #
    # F_Ω_p
    F_Ω_p .+= -R0T * (mul3(HT_p1,HT_p2,HT_p3,pdot) + 2/Δt*HT)

    # --- F_χ --- #
    # F_χ_V
    F_χ_V .+= F_χ_Vdot * 2/Δt
    # F_χ_Ω
    F_χ_Ω .+= F_χ_Ωdot * 2/Δt
    # F_χ_χ
    @inbounds for i in diagind(F_χ_χ)
        F_χ_χ[i] += 2/Δt
    end

end

# Inserts the element's Jacobian contributions into the global Jacobian matrix using precomputed maps for the indices
function insert_element_jacobians!(jacobian,F_u1_p,F_u2_p,F_u1_F,F_u2_F,F_u1_V,F_u2_V,F_u1_Ω,F_u2_Ω,F_u1_χ,F_u2_χ,F_u1_δ,F_u2_δ,F_p1_p,F_p2_p,F_p1_F,F_p2_F,F_p1_M,F_p2_M,F_p1_V,F_p2_V,F_p1_Ω,F_p2_Ω,F_p1_χ,F_p2_χ,F_p1_δ,F_p2_δ,F_F1_u,F_F2_u,F_F1_p,F_F2_p,F_F1_F,F_F2_F,F_F1_M,F_F2_M,F_M1_p,F_M2_p,F_M1_F,F_M2_F,F_M1_M,F_M2_M,F_V_u,F_V_p,F_V_V,F_Ω_p,F_Ω_Ω,F_χ_p,F_χ_V,F_χ_Ω,F_χ_χ,forceScaling,invForceScaling,hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat,jacobianMaps)

    # Unpack
    @unpack Fu1_p,Fu1_F,Fu1_V,Fu1_Ω,Fu1_χ,Fu1_δ,Fp1_p,Fp1_F,Fp1_M,Fp1_V,Fp1_Ω,Fp1_χ,Fp1_δ,FF1_u,FF1_p,FF1_F,FF1_M,FM1_p,FM1_F,FM1_M,Fu2_p,Fu2_F,Fu2_V,Fu2_Ω,Fu2_χ,Fu2_δ,Fp2_p,Fp2_F,Fp2_M,Fp2_V,Fp2_Ω,Fp2_χ,Fp2_δ,FF2_u,FF2_p,FF2_F,FF2_M,FM2_p,FM2_F,FM2_M,FV_u,FV_p,FV_V,FΩ_p,FΩ_Ω,Fχ_p,Fχ_V,Fχ_Ω,Fχ_χ = jacobianMaps
    @unpack nzval = jacobian

    ## Node 1 - equilibrium
    # --------------------------------------------------------------------------
    # --- F_u1 --- #
    set_block!(nzval, Fu1_p, F_u1_p * invForceScaling)
    set_block!(nzval, Fu1_F, F_u1_F)
    set_block!(nzval, Fu1_V, F_u1_V * invForceScaling)
    set_block!(nzval, Fu1_Ω, F_u1_Ω * invForceScaling)
    set_block!(nzval, Fu1_χ, F_u1_χ * invForceScaling)
    add_block!(nzval, Fu1_δ, F_u1_δ * invForceScaling)
    # --- F_p1 --- #
    set_block!(nzval, Fp1_p, notHingedNode1Mat * F_p1_p * invForceScaling)
    set_block!(nzval, Fp1_F, notHingedNode1Mat * F_p1_F)
    set_block!(nzval, Fp1_M, notHingedNode1Mat * F_p1_M)
    set_block!(nzval, Fp1_V, notHingedNode1Mat * F_p1_V * invForceScaling)
    set_block!(nzval, Fp1_Ω, notHingedNode1Mat * F_p1_Ω * invForceScaling)
    set_block!(nzval, Fp1_χ, F_p1_χ * invForceScaling)
    add_block!(nzval, Fp1_δ, F_p1_δ * invForceScaling)

    ## Node 1 — compatibility
    # --------------------------------------------------------------------------
    # --- F_F1 --- #
    set_block!(nzval, FF1_u, F_F1_u)
    set_block!(nzval, FF1_p, F_F1_p)
    set_block!(nzval, FF1_F, F_F1_F * forceScaling)
    set_block!(nzval, FF1_M, F_F1_M * forceScaling)
    # --- F_M1 --- #
    set_block!(nzval, FM1_p, notHingedNode1Mat * F_M1_p + hingedNode1Mat * F_p1_p * invForceScaling)
    set_block!(nzval, FM1_F, notHingedNode1Mat * F_M1_F * forceScaling + hingedNode1Mat * F_p1_F)
    set_block!(nzval, FM1_M, notHingedNode1Mat * F_M1_M * forceScaling + hingedNode1Mat * F_p1_M)

    ## Node 2 — equilibrium
    # --------------------------------------------------------------------------
    # --- F_u2 --- #
    set_block!(nzval, Fu2_p, F_u2_p * invForceScaling)
    set_block!(nzval, Fu2_F, F_u2_F)
    set_block!(nzval, Fu2_V, F_u2_V * invForceScaling)
    set_block!(nzval, Fu2_Ω, F_u2_Ω * invForceScaling)
    set_block!(nzval, Fu2_χ, F_u2_χ * invForceScaling)
    add_block!(nzval, Fu2_δ, F_u2_δ * invForceScaling)
    # --- F_p2 --- #
    set_block!(nzval, Fp2_p, F_p2_p * invForceScaling)
    set_block!(nzval, Fp2_F, F_p2_F)
    set_block!(nzval, Fp2_M, F_p2_M)
    set_block!(nzval, Fp2_V, F_p2_V * invForceScaling)
    set_block!(nzval, Fp2_Ω, F_p2_Ω * invForceScaling)
    set_block!(nzval, Fp2_χ, F_p2_χ * invForceScaling)
    add_block!(nzval, Fp2_δ, F_p2_δ * invForceScaling)

    ## Node 2 — compatibility
    # --------------------------------------------------------------------------
    # --- F_F2 --- #
    set_block!(nzval, FF2_u, F_F2_u)
    set_block!(nzval, FF2_p, F_F2_p)
    set_block!(nzval, FF2_F, F_F2_F * forceScaling)
    set_block!(nzval, FF2_M, F_F2_M * forceScaling)
    # --- F_M2 --- #
    set_block!(nzval, FM2_p, notHingedNode2Mat * F_M2_p)
    set_block!(nzval, FM2_F, notHingedNode2Mat * F_M2_F * forceScaling)
    set_block!(nzval, FM2_M, notHingedNode2Mat * F_M2_M * forceScaling)

    ## Midpoint - generalized velocities
    # --------------------------------------------------------------------------
    # --- F_V --- #
    set_block!(nzval, FV_u, F_V_u)
    set_block!(nzval, FV_p, F_V_p)
    set_block!(nzval, FV_V, F_V_V)
    # --- F_Ω --- #
    set_block!(nzval, FΩ_p, F_Ω_p)
    set_block!(nzval, FΩ_Ω, F_Ω_Ω)
    # --- F_χ --- #
    set_block!(nzval, Fχ_p, F_χ_p)
    set_block!(nzval, Fχ_V, F_χ_V)
    set_block!(nzval, Fχ_Ω, F_χ_Ω)
    set_block!(nzval, Fχ_χ, F_χ_χ)

end

# Computes the contributions from the current element to the inertia matrix
function element_inertia!(inertia,model,element,timeNow,Δt,isDynamicProblem)

    @unpack invForceScaling = model.assemblyData
    @unpack inertiaMaps = element
    @unpack Δℓ,R0,R0T = element.geometry
    @unpack I11,I12,I21,I22 = element.sectional
    @unpack P,H = element.state.midpoint
    @unpack RR0,HT,R_p1,R_p2,R_p3 = element.stateCache.rotationMatrices
    @unpack f1χ_Vdot,f2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Vdot,m2χ_Vdot,m1χ_Ωdot,m2χ_Ωdot = element.stateCache.aeroLoadsJacobians
    @unpack F_u1_pdot,F_u2_pdot,F_u1_Vdot,F_u2_Vdot,F_u1_Ωdot,F_u2_Ωdot,F_p1_pdot,F_p2_pdot,F_p1_Vdot,F_p2_Vdot,F_p1_Ωdot,F_p2_Ωdot,F_V_udot,F_Ω_pdot,F_χ_Vdot,F_χ_Ωdot,F_χ_χdot = element.stateCache.jacobians

    # Derivatives of the rotation variables w.r.t. the rotation parameters (time derivative of rotation tensor not required)
    element_rotation_derivatives!(element,Δt,false)

    # Compute aerodynamic contributions
    aero_derivatives_states_rates!(model,element,timeNow,Δt,isDynamicProblem)

    # Compute resultants
    element_inertia_resultants!(F_u1_pdot,F_u2_pdot,F_u1_Vdot,F_u2_Vdot,F_u1_Ωdot,F_u2_Ωdot,F_p1_pdot,F_p2_pdot,F_p1_Vdot,F_p2_Vdot,F_p1_Ωdot,F_p2_Ωdot,F_V_udot,F_Ω_pdot,Δℓ,R0,R0T,I11,I12,I21,I22,P,H,RR0,HT,R_p1,R_p2,R_p3,f1χ_Vdot,f2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Vdot,m2χ_Vdot,m1χ_Ωdot,m2χ_Ωdot)

    # Insert resultants into the inertia matrix
    insert_element_inertias!(inertia,F_u1_pdot,F_u2_pdot,F_u1_Vdot,F_u2_Vdot,F_u1_Ωdot,F_u2_Ωdot,F_p1_pdot,F_p2_pdot,F_p1_Vdot,F_p2_Vdot,F_p1_Ωdot,F_p2_Ωdot,F_V_udot,F_Ω_pdot,F_χ_Vdot,F_χ_Ωdot,F_χ_χdot,invForceScaling,inertiaMaps)

    # @logmsg LogLevel(show_debug() ? 0 : -1) "ELEMENT $(element.assemblyData.globalID) INERTIAS"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_u1_pdot: $F_u1_pdot"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_u1_Vdot: $F_u1_Vdot"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_u1_Ωdot: $F_u1_Ωdot"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p1_pdot: $F_p1_pdot"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p1_Vdot: $F_p1_Vdot"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p1_Ωdot: $F_p1_Ωdot"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_V_udot: $F_V_udot"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_Ω_pdot: $F_Ω_pdot"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_χ_Vdot: $F_χ_Vdot"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_χ_Ωdot: $F_χ_Ωdot"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_χ_χdot: $F_χ_χdot"
    # @logmsg LogLevel(show_debug() ? 0 : -1) (norm(inertia))

end

# Computes the derivatives of the aerodynamic loads and aerodynamic states w.r.t elemental states' rates
function aero_derivatives_states_rates!(model,element,timeNow,Δt,isDynamicProblem)

    # Skip if there are no aero loads, or angle of attack is undefined
    if isnothing(element.aero.surface) || isnan(element.aero.stateCache.attachedVars.steadyKin.α)
        @unpack f1χ_Vdot,f2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Vdot,m2χ_Vdot,m1χ_Ωdot,m2χ_Ωdot = element.stateCache.aeroLoadsJacobians
        @unpack F_χ_Vdot,F_χ_Ωdot = element.stateCache.jacobians
        zero_arrays!(f1χ_Vdot,f2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Vdot,m2χ_Vdot,m1χ_Ωdot,m2χ_Ωdot,F_χ_Vdot,F_χ_Ωdot)
        return
    end

    # Unpack
    @unpack derivationMethod = element.aero.surface.solvers
    @unpack Vdot,Ωdot = element.state.rates

    # Set input states' rates for aero loads wrapper function
    rates = element.aero.stateCache.statesRatesBuffer
    copyto!(view(rates,1:3), Vdot)
    copyto!(view(rates,4:6), Ωdot)

    # Jacobians of aerodynamic loads and state matrices w.r.t. states' rates
    aero_jacobians_states_rates!(rates,model,element,timeNow,Δt,isDynamicProblem,derivationMethod)

end

# Computes the aerodynamic Jacobians of states' rates with Automatic Differentiation
function aero_jacobians_states_rates!(rates,model,element,timeNow,Δt,isDynamicProblem,derivationMethod::AD)

    # Derivatives of aerodynamic loads and state matrices w.r.t. states' rates
    J = ForwardDiff.jacobian(x -> wrapper_aero_loads_and_matrices_from_states_rates!(x,model,element,timeNow,Δt,isDynamicProblem), SVector(rates))

    # Set derivatives
    set_aero_jacobians_of_states_rates!(J,element)

end

# Computes the aerodynamic Jacobians of states' rates with Finite Differences
function aero_jacobians_states_rates!(rates,model,element,timeNow,Δt,isDynamicProblem,derivationMethod::FD)

    # Derivatives of aerodynamic loads and state matrices w.r.t. states' rates
    J = first(FiniteDifferences.jacobian(derivationMethod.method, x -> wrapper_aero_loads_and_matrices_from_states_rates!(x,model,element,timeNow,Δt,isDynamicProblem), rates))

    # Set derivatives
    set_aero_jacobians_of_states_rates!(J,element)

end

# Wrapper function for aerodynamic loads and state matrices using elemental states' rates as inputs
function wrapper_aero_loads_and_matrices_from_states_rates!(rates,model,element,timeNow,Δt,isDynamicProblem)

    # Unpack states
    @unpack p,V,Ω,χ = element.state.midpoint
    @unpack δ = element.aero.surface.state
    Vdot = @view rates[1:3]
    Ωdot = @view rates[4:6]

    # Set state cache for derivatives calculations
    stateCacheDiff = get_diff_state_cache!(element.aero,eltype(rates))

    # Unpack resultants and arrays
    @unpack f1χ,f2χ,m1χ,m2χ = stateCacheDiff.loadsResultants
    @unpack Ac,Ai,An,Af,Ag,A,B = stateCacheDiff.stateMatrices

    # Nodal resutants and state arrays
    f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,B = aero_core!(f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,A,B,model,element,stateCacheDiff,p,V,Ω,χ,δ,Vdot,Ωdot,timeNow,Δt,isDynamicProblem)

    # Set outputs array
    outBuffer = stateCacheDiff.outBuffer
    write_aero_jacobians_outputs!(outBuffer, f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,B)

    return outBuffer

end

# Computes the aerodynamic derivatives of the element's states' rates
function set_aero_jacobians_of_states_rates!(J,element)

    # Unpack element data
    @unpack χ = element.state.midpoint
    @unpack f1χ_Vdot,f2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Vdot,m2χ_Vdot,m1χ_Ωdot,m2χ_Ωdot = element.stateCache.aeroLoadsJacobians
    @unpack F_χ_Vdot,F_χ_Ωdot = element.stateCache.jacobians
    @unpack ppSolver = element.aero.surface.solvers
    @unpack Ac,Ai,An,Af,Ag = element.aero.stateCache.stateMatrices
    @unpack Ac_Vdot,Ai_Vdot,An_Vdot,Af_Vdot,Ag_Vdot,Ac_Ωdot,Ai_Ωdot,An_Ωdot,Af_Ωdot,Ag_Ωdot,A_Vdot,A_Ωdot,B_Vdot,B_Ωdot = element.aero.matricesDerivatives

    # Derivatives of aerodynamic nodal loads resultants w.r.t. Vdot
    f1χ_Vdot .= @view J[1:3, 1:3]
    f2χ_Vdot .= @view J[4:6, 1:3]
    m1χ_Vdot .= @view J[7:9, 1:3]
    m2χ_Vdot .= @view J[10:12, 1:3]

    # Derivatives of aerodynamic nodal loads resultants w.r.t. Ωdot
    f1χ_Ωdot .= @view J[1:3, 4:6]
    f2χ_Ωdot .= @view J[4:6, 4:6]
    m1χ_Ωdot .= @view J[7:9, 4:6]
    m2χ_Ωdot .= @view J[10:12, 4:6]

    # Derivatives of aerodynamic equations w.r.t. Vdot and Ωdot
    if !isempty(χ)

        # Set lengths of arrays
        Nc = length(Ac)
        Ni = length(Ai)
        Nn = length(An)
        Nf = length(Af)
        Ng = length(Ag)

        # Derivatives of state arrays
        Ac_Vdot .= @view J[13 : 12+Nc, 1:3]
        Ac_Ωdot .= @view J[13 : 12+Nc, 4:6]
        Ai_Vdot .= @view J[13+Nc : 12+Nc+Ni, 1:3]
        Ai_Ωdot .= @view J[13+Nc : 12+Nc+Ni, 4:6]
        An_Vdot .= @view J[13+Nc+Ni : 12+Nc+Ni+Nn, 1:3]
        An_Ωdot .= @view J[13+Nc+Ni : 12+Nc+Ni+Nn, 4:6]
        Af_Vdot .= @view J[13+Nc+Ni+Nn : 12+Nc+Ni+Nn+Nf, 1:3]
        Af_Ωdot .= @view J[13+Nc+Ni+Nn : 12+Nc+Ni+Nn+Nf, 4:6]
        Ag_Vdot .= @view J[13+Nc+Ni+Nn+Nf : 12+Nc+Ni+Nn+Nf+Ng, 1:3]
        Ag_Ωdot .= @view J[13+Nc+Ni+Nn+Nf : 12+Nc+Ni+Nn+Nf+Ng, 4:6]

        # Set A_Vdot and A_Ωdot
        @inbounds for k=1:3
            @views begin
                i = 1
                # Set _c block
                if ppSolver isa Inflow
                    n = isqrt(Nc)
                    A_Vdot[i:i+n-1, i:i+n-1, k] .= reshape(Ac_Vdot[:,k],(n,n,1))
                    A_Ωdot[i:i+n-1, i:i+n-1, k] .= reshape(Ac_Ωdot[:,k],(n,n,1))
                    i += n
                else
                    i = assemble_A_diagonal_block!(A_Vdot,A_Ωdot,Ac_Vdot[:,k],Ac_Ωdot[:,k],i,k)
                end
                # Set _i, _n, _f and _g diagonal blocks
                i = assemble_A_diagonal_block!(A_Vdot,A_Ωdot,Ai_Vdot[:,k],Ai_Ωdot[:,k],i,k)
                i = assemble_A_diagonal_block!(A_Vdot,A_Ωdot,An_Vdot[:,k],An_Ωdot[:,k],i,k)
                i = assemble_A_diagonal_block!(A_Vdot,A_Ωdot,Af_Vdot[:,k],Af_Ωdot[:,k],i,k)
                i = assemble_A_diagonal_block!(A_Vdot,A_Ωdot,Ag_Vdot[:,k],Ag_Ωdot[:,k],i,k)
            end
        end

        # Set B_Vdot and B_Ωdot
        B_Vdot .= @view J[end-length(χ)+1:end, 1:3]
        B_Ωdot .= @view J[end-length(χ)+1:end, 4:6]

        # Set Jacobians
        @inbounds for i=1:3
            F_χ_Vdot[:,i] .= -(A_Vdot[:,:,i] * χ + B_Vdot[:,i])
            F_χ_Ωdot[:,i] .= -(A_Ωdot[:,:,i] * χ + B_Ωdot[:,i])
        end

    end

end

# Assembles a diagonal block in A_Vdot and A_Ωdot
@inline function assemble_A_diagonal_block!(A_V,A_Ω,v1,v2,i,k)
    n = length(v1)
    @inbounds for j = 1:n
        A_V[i+j-1, i+j-1, k] = v1[j]
        A_Ω[i+j-1, i+j-1, k] = v2[j]
    end
    return i + n
end

# Computes the structural and aerodynamic contributions from the current element to the inertia matrix
function element_inertia_resultants!(F_u1_pdot,F_u2_pdot,F_u1_Vdot,F_u2_Vdot,F_u1_Ωdot,F_u2_Ωdot,F_p1_pdot,F_p2_pdot,F_p1_Vdot,F_p2_Vdot,F_p1_Ωdot,F_p2_Ωdot,F_V_udot,F_Ω_pdot,Δℓ,R0,R0T,I11,I12,I21,I22,P,H,RR0,HT,R_p1,R_p2,R_p3,f1χ_Vdot,f2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Vdot,m2χ_Vdot,m1χ_Ωdot,m2χ_Ωdot)

    # --- F_u --- #
    # F_u_pdot
    tmp = Δℓ/2 * mul3(R_p1,R_p2,R_p3,SVector(R0*P))
    F_u1_pdot .= tmp
    F_u2_pdot .= tmp
    # F_u_Vdot
    tmp = Δℓ/2 * SMatrix(RR0 * I11)
    F_u1_Vdot .= tmp - f1χ_Vdot
    F_u2_Vdot .= tmp - f2χ_Vdot
    # F_u_Ωdot
    tmp = Δℓ/2 * SMatrix(RR0 * I12)
    F_u1_Ωdot .= tmp - f1χ_Ωdot
    F_u2_Ωdot .= tmp - f2χ_Ωdot

    # --- F_p --- #
    # F_p_pdot
    tmp = Δℓ/2 * mul3(R_p1,R_p2,R_p3,SVector(R0*H))
    F_p1_pdot .= tmp
    F_p2_pdot .= tmp
    # F_p_Vdot
    tmp = Δℓ/2 * SMatrix(RR0 * I21)
    F_p1_Vdot .= tmp - m1χ_Vdot
    F_p2_Vdot .= tmp - m2χ_Vdot
    # F_p_Ωdot
    tmp = Δℓ/2 * SMatrix(RR0 * I22)
    F_p1_Ωdot .= tmp - m1χ_Ωdot
    F_p2_Ωdot .= tmp - m2χ_Ωdot

    # --- F_V --- #
    # F_V_udot (already set as constant)
    # F_V_udot .= - I3

    # --- F_Ω --- #
    # F_Ω_pdot
    F_Ω_pdot .= - SMatrix(R0T * HT)

end

# Inserts the element's contributions to the inertia matrix into the global inertia matrix
function insert_element_inertias!(inertia,F_u1_pdot,F_u2_pdot,F_u1_Vdot,F_u2_Vdot,F_u1_Ωdot,F_u2_Ωdot,F_p1_pdot,F_p2_pdot,F_p1_Vdot,F_p2_Vdot,F_p1_Ωdot,F_p2_Ωdot,F_V_udot,F_Ω_pdot,F_χ_Vdot,F_χ_Ωdot,F_χ_χdot,invForceScaling,inertiaMaps)

    @unpack Fu1_pdot,Fu1_Vdot,Fu1_Ωdot,Fp1_pdot,Fp1_Vdot,Fp1_Ωdot,Fu2_pdot,Fu2_Vdot,Fu2_Ωdot,Fp2_pdot,Fp2_Vdot,Fp2_Ωdot,FV_udot,FΩ_pdot,Fχ_Vdot,Fχ_Ωdot,Fχ_χdot = inertiaMaps
    @unpack nzval = inertia

    ## Node 1
    # --------------------------------------------------------------------------
    # --- F_u1 --- #
    set_block!(nzval, Fu1_pdot, F_u1_pdot * invForceScaling)
    set_block!(nzval, Fu1_Vdot, F_u1_Vdot * invForceScaling)
    set_block!(nzval, Fu1_Ωdot, F_u1_Ωdot * invForceScaling)
    # --- F_p1 --- #
    set_block!(nzval, Fp1_pdot, F_p1_pdot * invForceScaling)
    set_block!(nzval, Fp1_Vdot, F_p1_Vdot * invForceScaling)
    set_block!(nzval, Fp1_Ωdot, F_p1_Ωdot * invForceScaling)

    ## Node 2
    # --------------------------------------------------------------------------
    # --- F_u2 --- #
    set_block!(nzval, Fu2_pdot, F_u2_pdot * invForceScaling)
    set_block!(nzval, Fu2_Vdot, F_u2_Vdot * invForceScaling)
    set_block!(nzval, Fu2_Ωdot, F_u2_Ωdot * invForceScaling)
    # --- F_p2 --- #
    set_block!(nzval, Fp2_pdot, F_p2_pdot * invForceScaling)
    set_block!(nzval, Fp2_Vdot, F_p2_Vdot * invForceScaling)
    set_block!(nzval, Fp2_Ωdot, F_p2_Ωdot * invForceScaling)

    ## Midpoint
    # --------------------------------------------------------------------------
    # --- F_V --- #
    set_block!(nzval, FV_udot, F_V_udot)
    # --- F_Ω --- #
    set_block!(nzval, FΩ_pdot, F_Ω_pdot)
    # --- F_χ --- #
    set_block!(nzval, Fχ_Vdot, F_χ_Vdot)
    set_block!(nzval, Fχ_Ωdot, F_χ_Ωdot)
    set_block!(nzval, Fχ_χdot, F_χ_χdot)

end

# Computes the nodal states
function special_node_states!(model,specialNode,x,σ)

    @unpack forceScaling = model.assemblyData
    @unpack BCs = specialNode
    @unpack uIsPrescribed,pIsPrescribed,trimIsPrescribed = specialNode.flags
    @unpack u,p,F,M = specialNode.state
    @unpack DOF_uF,DOF_pM,DOF_T = specialNode.assemblyData

    # Check if node is not BC'ed - all nodal states are generalized displacements, externally applied forces are null
    if isempty(BCs)
        u .= x[DOF_uF]
        p .= x[DOF_pM]
        F .= 0.0
        M .= 0.0
        return
    end

    # Reset loads to zero
    zero_arrays!(F,M)

    # Loop boundary conditions
    for BC in BCs

        @unpack isLoad,isFollower,isTrim = BC.deducedData
        @unpack currentValue = BC.state
        @unpack F̂,M̂,Ftrim,Mtrim = BC.stateCache

        ## Generalized displacements
        # ----------------------------------------------------------------------
        # Loop directions and set displacements and rotations, resolved in basis A, as either nodal states or prescribed values
        @inbounds for i=1:3
            if uIsPrescribed[i] && !isLoad[i]
                # u[i] is a prescribed value
                u[i] = σ * currentValue[i]
            elseif !uIsPrescribed[i] && isLoad[i]
                # u[i] is a nodal state
                u[i] = x[DOF_uF[i]]
            end
            if pIsPrescribed[i] && !isLoad[i+3]
                # p[i] is a prescribed value
                p[i] = σ * currentValue[i+3]
            elseif !pIsPrescribed[i] && isLoad[i+3]
                # p[i] is a nodal state 
                p[i] = x[DOF_pM[i]]
            end
        end

        ## Prescribed values of generalized forces, resolved in basis A
        # ----------------------------------------------------------------------
        # Reset prescribed (hatted) and trim nodal values of current BC
        zero_arrays!(F̂,M̂,Ftrim,Mtrim)

        # Nodal rotation tensor from basis b to basis B
        R = rotation_tensor_WM(p)
        
        # Loop orthogonal directions
        @inbounds for i=1:3
            # Add non-trim forces
            if isLoad[i] && !isTrim[i]
                # Add dead forces (dead values are constant when resolved in basis A)
                F̂ .+= σ * currentValue[i] * !isFollower[i] * I3[:,i]
                # Add follower forces (follower values are constant when resolved in basis B - transform to basis A)
                F̂ .+= σ * currentValue[i] * isFollower[i] * R[:,i]
            end
            # Add trim forces
            if isTrim[i]
                # Add dead trim forces
                Ftrim .+= x[DOF_T[i]] * !isFollower[i] * I3[:,i] * forceScaling
                # Add follower trim forces
                Ftrim .+= x[DOF_T[i]] * isFollower[i] * R[:,i] * forceScaling
            end
            # Add non-trim moments
            if isLoad[i+3] && !isTrim[i+3]
                # Add dead moments (dead values are constant when resolved in basis A)
                M̂ .+= σ * currentValue[i+3] * !isFollower[i+3] * I3[:,i]
                # Add follower moments (follower values are constant when resolved in basis B - transform to basis A)
                M̂ .+= σ * currentValue[i+3] * isFollower[i+3] * R[:,i]
            end
            # Add trim moments
            if isTrim[i+3]
                # Add dead trim moments
                Mtrim .+= x[DOF_T[i+3]] * !isFollower[i+3] * I3[:,i] * forceScaling
                # Add follower trim moments
                Mtrim .+= x[DOF_T[i+3]] * isFollower[i+3] * R[:,i] * forceScaling
            end
        end
    
        ## Generalized forces ------------------------------------------------------------------------
        # Loop directions and set forces and moments, resolved in basis A, as either nodal states or prescribed values
        @inbounds for i=1:3
            if isLoad[i] && !isTrim[i]
                # F[i] is a prescribed value (plus the component of trim follower forces, which rotate into this direction)
                F[i] += F̂[i] + Ftrim[i]
            elseif isTrim[i]
                # F[i] is a trim value
                F[i] = Ftrim[i]
            elseif !trimIsPrescribed[i]
                # F[i] is a nodal state
                F[i] = x[DOF_uF[i]] * forceScaling
            end
            if isLoad[i+3] && !isTrim[i+3]
                # M[i] is a prescribed value (plus the component of trim follower moments, which rotate into this direction)
                M[i] += M̂[i] + Mtrim[i]
            elseif isTrim[i+3] 
                # M[i] is a trim value
                M[i] = Mtrim[i]
            elseif !trimIsPrescribed[i+3]
                # M[i] is a nodal state
                M[i] = x[DOF_pM[i]] * forceScaling
            end
        end
    end

end

# Adds the contributions of attached springs to the special node's loads
function spring_loads!(model,specialNode)

    @unpack springs = specialNode
    @unpack doubleAttachmentSpringsData = specialNode.assemblyData
    @unpack u,p,F,M = specialNode.state

    # Initialize counter for double attachment springs
    c = 0

    # Loop springs
    for spring in springs
        @unpack hasDoubleAttachment = spring.assemblyData
        @unpack FsFun,MsFun = spring.assemblyConstitutive
        @unpack Δu,Δp,Fs,Ms = spring.state
        # Generalized spring displacements, resolved in basis A
        if hasDoubleAttachment # Double attachment springs
            # Update counter
            c += 1
            # Unpack double attachment spring data for this spring
            @unpack otherNodeSpecialID,uOtherNode,pOtherNode,otherNodeDOF_uF,otherNodeDOF_pM = doubleAttachmentSpringsData[c]
            # Node of other attachment point
            otherNode = model.assemblyData.specialNodes[otherNodeSpecialID]
            # Other node's states
            uOtherNode .= otherNode.state.u
            pOtherNode .= otherNode.state.p
            # Generalized spring displacements
            Δu .= u - uOtherNode
            Δp .= rotation_between_WM(pOtherNode,p)
        else # Single attachment springs
            Δu .= u
            Δp .= scaled_rotation_parameters(p)
        end
        # Spring loads, resolved in basis A
        Fs .= - FsFun(Δu)
        Ms .= - MsFun(Δp)
        # Add spring loads to nodal loads
        F .+= Fs
        M .+= Ms
    end

end

# Computes the contributions from the current node to the residual array
function special_node_residual!(residual,model,specialNode)

    @unpack invForceScaling = model.assemblyData
    @unpack springs = specialNode
    @unpack u,p,F,M = specialNode.state
    @unpack ζonElements,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM = specialNode.assemblyData

    # @logmsg LogLevel(show_debug() ? 0 : -1) "SPECIAL NODE $(specialNode.assemblyData.globalID) RESIDUALS"

    # Add spring loads, if applicable
    isempty(springs) || spring_loads!(model,specialNode)

    # Insert the contributions from the current special node into the residual array
    insert_special_node_residual!(residual,u,p,F,M,ζonElements,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,invForceScaling)

    # @logmsg LogLevel(show_debug() ? 0 : -1) residual[eqs_Fu]
    # @logmsg LogLevel(show_debug() ? 0 : -1) residual[eqs_Fp]
    # for e in eachindex(ζonElements)
    #     @logmsg LogLevel(show_debug() ? 0 : -1) residual[eqs_FF[e]]
    #     @logmsg LogLevel(show_debug() ? 0 : -1) residual[eqs_FM[e]]
    # end

end

# Inserts the contributions from the current special node into the residual array
function insert_special_node_residual!(residual,u,p,F,M,ζonElements,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,invForceScaling)

    ## Equilibrium equations
    # ----------------------------------------------------------------------
    # Add to node's equilibrium and compatibility equations' residuals
    residual[eqs_Fu] .-= F * invForceScaling
    residual[eqs_Fp] .-= M * invForceScaling

    ## Compatibility equations
    # ----------------------------------------------------------------------
    # Loop connected elements of the node
    for e in eachindex(ζonElements)
        # Add to node's compatibility equations' residuals
        residual[eqs_FF[e]] .+= ζonElements[e] * u
        residual[eqs_FM[e]] .+= ζonElements[e] * p
    end
    
end

# Computes the contributions from the current node to the Jacobian matrix
function special_node_jacobian!(jacobian,model,specialNode,x,σ)

    @unpack invForceScaling = model.assemblyData
    @unpack BCs,springs,jacobianMaps = specialNode
    @unpack uIsPrescribed,pIsPrescribed,allBCflags,uNotPrescribedMat,pNotPrescribedMat = specialNode.flags
    @unpack p,R,F_p,M_p = specialNode.state
    @unpack ζonElements,doubleAttachmentSpringsData,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,DOF_uF,DOF_pM,DOF_T = specialNode.assemblyData

    # @logmsg LogLevel(show_debug() ? 0 : -1) "SPECIAL NODE $(specialNode.assemblyData.globalID) JACOBIANS"

    # Follower loads' derivatives w.r.t. the rotation parameters (contributions to Jacobian matrix)
    special_node_follower_loads_derivatives_rotation_parameters!(jacobian,model,specialNode,x,σ)

    # Add contributions from the current BC (if any) on the current node
    for (isLoad, isTrim, isFollower) in (isempty(BCs) ? _emptyBCflags : allBCflags)
        insert_special_node_jacobian!(jacobian,invForceScaling,F_p,M_p,R,ζonElements,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,DOF_uF,DOF_pM,DOF_T,uIsPrescribed,pIsPrescribed,isLoad,isTrim,isFollower)
    end

    # if !isempty(springs)
    #     @logmsg LogLevel(show_debug() ? 0 : -1) "SPECIAL NODE $(specialNode.assemblyData.globalID) SPRING JACOBIANS"
    # end

    # Add contributions from the springs attached to the current node
    spring_loads_jacobians!(model,springs,jacobian,jacobianMaps,invForceScaling,doubleAttachmentSpringsData,p,uNotPrescribedMat,pNotPrescribedMat)

end

# Computes the contributions from the nodal follower loads to the Jacobian matrix
function special_node_follower_loads_derivatives_rotation_parameters!(jacobian,model,specialNode,x,σ)

    @unpack forceScaling = model.assemblyData
    @unpack BCs = specialNode
    @unpack hasFollowers = specialNode.flags
    @unpack p,R,F_p,M_p = specialNode.state
    @unpack DOF_T = specialNode.assemblyData

    # Skip if there are no follower loads
    hasFollowers || return

    # Reset to zero
    zero_arrays!(F_p,M_p)

    # Nodal rotation tensor from basis b to basis B (gives the directions of the follower loads)
    R .= rotation_tensor_WM(p)

    # Get nodal rotation tensor's derivatives w.r.t. the extended rotation parameters
    R_p1,R_p2,R_p3 = ∂R_∂p(p)

    # Loop boundary conditions
    for BC in BCs

        @unpack isLoad,isFollower,isTrim = BC.deducedData
        @unpack currentValue = BC.state
        @unpack Ff,Mf = BC.stateCache

        # Skip if there are no follower loads
        any(isFollower) || continue

        # Nodal follower load components, resolved in basis A
        @inbounds for i=1:3
            # Forces
            if isTrim[i] && isFollower[i]
                # Trim variable 
                Ff[i] = x[DOF_T[i]] * forceScaling
            elseif isFollower[i]
                # Not a trim variable
                Ff[i] = σ * currentValue[i]
            end
            # Moments
            if isTrim[i+3] && isFollower[i+3]
                # Trim variable 
                Mf[i] = x[DOF_T[i+3]] * forceScaling
            elseif isFollower[i+3]
                # Not a trim variable
                Mf[i] = σ * currentValue[i+3]
            end
        end

        # Add to derivatives of the follower loads w.r.t. the extended rotation parameters 
        F_p .+= mul3(R_p1,R_p2,R_p3,Ff)
        M_p .+= mul3(R_p1,R_p2,R_p3,Mf)

    end

end

# Updates the Jacobian matrix with the contributions from the current BC (if any) on the current node 
function insert_special_node_jacobian!(jacobian,invForceScaling,F_p,M_p,R,ζonElements,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,DOF_uF,DOF_pM,DOF_T,uIsPrescribed,pIsPrescribed,isLoad,isTrim,isFollower)

    ## Equilibrium equations
    # ----------------------------------------------------------------------
    # Outer loop over node's DOFs
    @inbounds for i=1:3
        # If F[i] is prescribed, u[i] is the nodal state
        if isLoad[i] && !uIsPrescribed[i]
            # Inner loop over node's DOFs
            @inbounds for j=1:3
                # If M[j] is prescribed, p[j] is the nodal state
                if isLoad[j+3] && !pIsPrescribed[j]
                    # Equilibrium Jacobian: d(-F)/d(p) = -∂F/∂p
                    jacobian[eqs_Fu[i],DOF_pM[j]] = -F_p[i,j] * invForceScaling
                end
            end
        # Otherwise u[i] is prescribed, and F[i] is the nodal state    
        else
            # Equilibrium Jacobian: d(-F)/d(F) = -1
            jacobian[eqs_Fu[i],DOF_uF[i]] = -1
        end
        # If F[i] is also a trim variable
        if isTrim[i]
            # Equilibrium Jacobian: d(-F)/d(T) = -1 along direction i for dead trim forces, or -R[:,i] for follower trim forces (which have components along all directions of basis A)
            @inbounds for j=1:3
                if isFollower[i] ? (isLoad[j] || isTrim[j]) : j == i
                    jacobian[eqs_Fu[j],DOF_T[i]] = isFollower[i] ? -R[j,i] : -1
                end
            end
        end
        # If M[i] is prescribed, p[i] is the nodal state
        if isLoad[i+3] && !pIsPrescribed[i]
            # Inner loop over node's DOFs
            @inbounds for j=1:3
                # If M[j] is prescribed, p[j] is the nodal state
                if isLoad[j+3] && !pIsPrescribed[j]
                    # Equilibrium Jacobian d(-M)/d(p) = -∂M/∂p
                    jacobian[eqs_Fp[i],DOF_pM[j]] = -M_p[i,j] * invForceScaling
                end
            end
        # If p[i] is prescribed, M[i] is the nodal state    
        else
            # Equilibrium Jacobian: d(-M)/d(M) = -1
            jacobian[eqs_Fp[i],DOF_pM[i]] = -1
        end
        # If M[i] is also a trim variable
        if isTrim[i+3]
            # Equilibrium Jacobian: d(-M)/d(T) = -1 along direction i for dead trim moments, or -R[:,i] for follower trim moments
            @inbounds for j=1:3
                if isFollower[i+3] ? (isLoad[j+3] || isTrim[j+3]) : j == i
                    jacobian[eqs_Fp[j],DOF_T[i+3]] = isFollower[i+3] ? -R[j,i] : -1
                end
            end
        end
    end

    ## Compatibility equations
    # ----------------------------------------------------------------------
    # Loop connected elements of the node
    for e in eachindex(ζonElements)
        # Loop over node's DOFs
        @inbounds for i=1:3
            # If F[i] is prescribed, u[i] is the nodal state
            if isLoad[i] && !uIsPrescribed[i]
                # Compatibility Jacobian: d(-/+u)/d(u) = -/+1 = ζonElements[e]
                jacobian[eqs_FF[e][i],DOF_uF[i]] = ζonElements[e]
            end
            # If M[i] is prescribed, p[i] is the nodal state
            if isLoad[i+3] && !pIsPrescribed[i]
                # Compatibility Jacobian: d(-/+p)/d(p) = -/+1 = ζonElements[e]
                jacobian[eqs_FM[e][i],DOF_pM[i]] = ζonElements[e]
            end
        end
    end

    # @logmsg LogLevel(show_debug() ? 0 : -1) Matrix(jacobian[[eqs_Fu; eqs_Fp], [DOF_uF; DOF_pM]])
    # for e in eachindex(ζonElements)
    #     @logmsg LogLevel(show_debug() ? 0 : -1) Matrix(jacobian[[eqs_FF[e]; eqs_FM[e]], [DOF_uF; DOF_pM]])
    # end
    # @logmsg LogLevel(show_debug() ? 0 : -1) (norm(jacobian))
    # if size(jacobian, 1) == size(jacobian, 2)
    #     @logmsg LogLevel(show_debug() ? 0 : -1) (det(jacobian))
    # end

end

# Adds the contributions of the spring loads to the Jacobian matrix
function spring_loads_jacobians!(model,springs,jacobian,jacobianMaps,invForceScaling,doubleAttachmentSpringsData,p,uNotPrescribedMat,pNotPrescribedMat)

    # Unpack node data
    @unpack Fu_uF,Fp_pM,Fu_uFother,Fp_pMother = jacobianMaps
    @unpack nzval = jacobian

    # Initialize counter for double attachment springs
    c = 0

    # Loop over springs
    for spring in springs

        # Unpack spring data
        @unpack hasDoubleAttachment = spring.assemblyData
        @unpack Fs_ΔuFun,Ms_ΔpFun = spring.assemblyConstitutive
        @unpack Δu,Δp = spring.state

        # Derivatives of spring loads w.r.t. deformations
        Fs_Δu = Fs_ΔuFun(Δu)
        Ms_Δp = Ms_ΔpFun(Δp)

        # Double attachment springs
        if hasDoubleAttachment
            # Update counter
            c += 1
            # Unpack double attachment spring data for this spring
            @unpack otherNodeSpecialID,otherNodeDOF_uF,otherNodeDOF_pM,pOtherNode = doubleAttachmentSpringsData[c]
            # Node of other attachment point
            otherNode = model.assemblyData.specialNodes[otherNodeSpecialID]
            # Add translational spring Jacobians (for DOFs that are not prescribed): Δu = u - uOtherNode → ∂(-Fs)/∂u = Fs_Δu and ∂(-Fs)/∂uOtherNode = -Fs_Δu
            add_block!(nzval, Fu_uF, Fs_Δu * invForceScaling * uNotPrescribedMat)
            add_block!(nzval, Fu_uFother[c], - Fs_Δu * invForceScaling * otherNode.flags.uNotPrescribedMat)
            # Derivative of spring rotation w.r.t rotation of each of its attachment nodes
            if isapprox(p, pOtherNode; atol=1e-10)
                Δp_p = I3
                Δp_pOtherNode = -I3
            else
                Δp_p = SMatrix{3,3}(ForwardDiff.jacobian(x -> rotation_between_WM(pOtherNode,x), p))
                Δp_pOtherNode = SMatrix{3,3}(ForwardDiff.jacobian(x -> rotation_between_WM(x,p), pOtherNode))
            end
            # Add rotational spring Jacobians (for DOFs that are not prescribed):  Δp = rotation_between_WM(p,pOtherNode) → ∂(-Ms)/∂p = ∂(-Ms)/∂Δp*∂Δp/∂p = Ms_Δp*Δp_p and ∂(-Ms)/∂pOtherNode = ∂(-Ms)/∂Δp*∂Δp/∂pOtherNode = Ms_Δp*Δp_pOtherNode
            add_block!(nzval, Fp_pM, Ms_Δp * Δp_p * invForceScaling * pNotPrescribedMat)
            add_block!(nzval, Fp_pMother[c], Ms_Δp * Δp_pOtherNode * invForceScaling * otherNode.flags.pNotPrescribedMat)
        # Single attachment springs
        else
            # Add translational spring Jacobians (for DOFs that are not prescribed): Δu = u → ∂(-Fs)/∂u = Fs_Δu
            add_block!(nzval, Fu_uF, Fs_Δu * invForceScaling * uNotPrescribedMat)
            # Add rotational spring Jacobian (for DOFs that are not prescribed):  Δp = scaled_rotation_parameters(p) → ∂(-Ms)/∂p = ∂(-Ms)/∂Δp*∂Δp/∂p = Ms_Δp*∂ps/∂p
            λ,_,_,pNorm,_ = rotation_parameter_scaling(p)
            Δp_p = ∂ps_∂p(λ,p,pNorm)
            add_block!(nzval, Fp_pM, Ms_Δp * Δp_p * invForceScaling * pNotPrescribedMat)
        end

    end

    # if !isempty(springs)
    #     @logmsg LogLevel(show_debug() ? 0 : -1) reshape(Vector(nzval[Fu_uF]), 3, 3)
    #     @logmsg LogLevel(show_debug() ? 0 : -1) reshape(Vector(nzval[Fp_pM]), 3, 3)
    #     if !isempty(Fu_uFother)
    #         @logmsg LogLevel(show_debug() ? 0 : -1) reshape(Vector(nzval[Fu_uFother[1]]), 3, 3)
    #         @logmsg LogLevel(show_debug() ? 0 : -1) reshape(Vector(nzval[Fp_pMother[1]]), 3, 3)
    #     end
    #     @logmsg LogLevel(show_debug() ? 0 : -1) (norm(jacobian))
    # end
    # if size(jacobian, 1) == size(jacobian, 2)
    #     @logmsg LogLevel(show_debug() ? 0 : -1) (det(jacobian))
    # end

end

# Computes the contributions from the hinge axis constraints to the residual array
function hac_residual!(residual,model,constraint,x)

    @unpack rotationIsFixed,initialHingeAxis,slaveDir,pHValuePrescribed = constraint
    @unpack eqs_Fp,eqs_λ = constraint.assemblyData
    @unpack λ,B,GᵀBᵀ,pM,pS = constraint.state

    # Update state
    hac_state!(model,constraint,x)
    
    # Insert residual contributions from the hinge axis constraint
    insert_hac_residual!(residual,rotationIsFixed,initialHingeAxis,slaveDir,pHValuePrescribed,eqs_Fp,eqs_λ,λ,GᵀBᵀ,pM,pS)

    # @logmsg LogLevel(show_debug() ? 0 : -1) "HAC STATES"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "λ: $λ"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "pM: $pM"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "pS: $pS"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "B: $B"
    
    # @logmsg LogLevel(show_debug() ? 0 : -1) "HAC RESIDUALS"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "C: $(residual[eqs_λ])"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "F_p: $(residual[eqs_Fp])"

end

# Updates the states of the hinge axis constraints
function hac_state!(model,constraint,x; updateOutputs=false)
    
    @unpack forceScaling = model.assemblyData
    @unpack rotationIsFixed,initialHingeAxis,slaveDir,pHValuePrescribed = constraint
    @unpack DOF_pM,DOF_pS,DOF_λ = constraint.assemblyData
    @unpack λ,B,Bᵀ,Gᵀ,GᵀBᵀ,pM,pS,balanceMoment = constraint.state

    # Lagrange multipliers
    λ .= x[DOF_λ]

    # Rotation parameters of master and slave elements
    pM .= x[DOF_pM]
    pS .= x[DOF_pS]

    # Mapping from virtual rotations to rotation parameters' variations
    Gᵀ .= GT(pM,pS)

    # Constraint matrix (Jacobian of the constraint equations) and its transpose
    B .= ∂C_∂p(pM,pS,initialHingeAxis,rotationIsFixed,pHValuePrescribed,slaveDir)
    Bᵀ .= B'

    # Transpose of the constraint matrix w.r.t. virtual rotations (constraint moments' matrix)
    GᵀBᵀ .= Gᵀ * Bᵀ

    # Update output states, if applicable
    if updateOutputs
        # Balance moment vector
        @views balanceMoment .= (GᵀBᵀ * λ)[4:6] * forceScaling
        # Hinge moment (magnitude of balance moment)
        hingeMoment = norm(balanceMoment)
        # Deformed hinge axis
        deformedHingeAxis = deformed_hinge_axis(pM,initialHingeAxis)
        # Hinge rotation parameters
        pH = hinge_rotation_parameters(pM,pS)
        pHValue = hinge_rotation_value(pH,deformedHingeAxis)
        # Hinge rotation angle
        ϕ = rotation_angle_limited(pH)
        @pack! constraint.state = deformedHingeAxis,pH,pHValue,ϕ,hingeMoment
    end

end

# Inserts the residual contributions from the hinge axis constraints into the residual vector
function insert_hac_residual!(residual,rotationIsFixed,initialHingeAxis,slaveDir,pHValuePrescribed,eqs_Fp,eqs_λ,λ,GᵀBᵀ,pM,pS)

    # Add to residual of nominal system
    residual[eqs_Fp] .+= GᵀBᵀ * λ

    # Set constraint equations' residual
    residual[eqs_λ] .= C(pM,pS,initialHingeAxis,rotationIsFixed,pHValuePrescribed,slaveDir)

end

# Computes the contributions from the hinge axis constraints to the Jacobian matrix
function hac_jacobian!(jacobian,model,constraint)

    @unpack invForceScaling = model.assemblyData
    @unpack rotationIsFixed,slaveDir,initialHingeAxis,pHValuePrescribed = constraint
    @unpack jacobianMaps = constraint
    @unpack λ,B,Bᵀ,Gᵀ,GᵀBᵀ,GᵀBᵀλ_p,pM,pS = constraint.state

    # Update constraint moments Jacobian w.r.t the master and slave elements' rotation parameters
    GᵀBᵀλ_p .= ∂GTBTλ_∂p(λ,Bᵀ,Gᵀ,pM,pS,initialHingeAxis,rotationIsFixed,pHValuePrescribed,slaveDir)

    # Insert the hinge axis constraints' contributions into the Jacobian matrix
    insert_hac_jacobian!(jacobian,jacobianMaps,B,GᵀBᵀ,GᵀBᵀλ_p)

    # @logmsg LogLevel(show_debug() ? 0 : -1) "HAC JACOBIANS"
    # @logmsg LogLevel(show_debug() ? 0 : -1) "GᵀBᵀλ_p: $GᵀBᵀλ_p"

end

# Inserts the hinge axis constraints' contributions into the Jacobian matrix
function insert_hac_jacobian!(jacobian,jacobianMaps,B,GᵀBᵀ,GᵀBᵀλ_p)

    @unpack nzval = jacobian
    @unpack Fp_pMS,Fp_λ,Fλ_pMS = jacobianMaps

    # --- F_p --- #
    add_block!(nzval, Fp_pMS, GᵀBᵀλ_p)
    set_block!(nzval, Fp_λ, GᵀBᵀ)

    # --- C --- #
    set_block!(nzval, Fλ_pMS, B)

end


# Updates the elemental and nodal states
function update_states!(problem)

    @unpack x = problem.systemArrays
    @unpack σ = problem.systemSolver.loadStepState
    @unpack model = problem
    timeNow = time_now(problem)
    timeIndex = time_index(problem)
    Δt = time_step(problem)
    isDynamicProblem = problem isa DynamicProblem

    update_states!(model,x,σ,timeIndex,timeNow,Δt,isDynamicProblem)

end
function update_states!(model,x,σ,timeIndex,timeNow,Δt,isDynamicProblem)

    @unpack hingeAxisConstraints = model
    @unpack elementGroups,specialNodes = model.assemblyData

    # Loop over elements
    foreach_element(elementGroups) do element

        ## Generalized velocities and accelerations of basis b at element's midpoint, resolved in basis A
        # ----------------------------------------------------------------------
        # Velocities
        element_velocities_basis_b!(model,element,σ)
        # Accelerations
        element_accelerations_basis_b!(model,element,σ)

        ## States and states' rates
        # ----------------------------------------------------------------------
        # States
        element_states!(model,element,x)
        # States' rates
        isDynamicProblem && element_states_rates!(element,Δt)
        
        ## Rotation variables
        # ----------------------------------------------------------------------
        element_rotation_variables!(element,isDynamicProblem)
        
        ## Distributed loads
        # ----------------------------------------------------------------------
        element_distributed_loads!(model,element,σ,timeIndex,timeNow,Δt,isDynamicProblem)

        ## Update nodal states of the element
        # ----------------------------------------------------------------------
        element_nodal_states!(element)

    end

    # Loop over hinge axis constraints
    for constraint in hingeAxisConstraints
        hac_state!(model,constraint,x; updateOutputs=true)
    end

    # Loop over special nodes
    for specialNode in specialNodes

        ## Update nodal states of the special node
        # ----------------------------------------------------------------------
        special_node_states!(model,specialNode,x,σ)

        ## Add contributions of attached springs to the special node's loads
        # ----------------------------------------------------------------------
        spring_loads!(model,specialNode)

    end

end

# Updates the nodal states of the element
function element_nodal_states!(element)

    @unpack Δℓ,k,R0,R0_n1,R0_n2,R0T_n1,R0T_n2 = element.geometry
    @unpack RR0,RR0T,HTinv = element.stateCache.rotationMatrices
    @unpack u,p,F,M,V,Ω,γ,κ,P,H = element.state.midpoint
    @unpack Pdot,Hdot = element.state.rates
    @unpack u1,u2,p1,p2,u1_b,u2_b,p1_b,p2_b,F1,F2,M1,M2,θ = element.state.nodal
    @unpack f1,f2,m1,m2 = element.stateCache.loadsResultants

    # Midpoint generalized displacements' derivatives w.r.t x1
    u_prime = RR0 * (a1 + γ) - R0 * a1
    p_prime = HTinv * SVector(R0 * κ)

    # Midpoint deformed curvature vector, resolved in basis B
    K = k + κ

    # Midpoint generalized forces' derivatives w.r.t x1
    F_prime = Pdot + cross(Ω,P) - cross(K,F)
    M_prime = Hdot + cross(Ω,H) + cross(V,P) - cross(a1+γ,F) - cross(K,M)

    # Nodal displacements and rotation parameters' vectors, resolved in basis A
    u1 .= u - Δℓ/2 * u_prime
    u2 .= u + Δℓ/2 * u_prime
    p1 .= p - Δℓ/2 * p_prime
    p2 .= p + Δℓ/2 * p_prime

    # Nodal displacements and rotation parameters' vectors, resolved in basis b
    u1_b .= R0T_n1 * u1
    u2_b .= R0T_n2 * u2
    p1_b .= R0T_n1 * p1
    p2_b .= R0T_n2 * p2

    # Nodal generalized internal forces, resolved in basis B
    F1 .= F - (Δℓ/2 * F_prime - RR0T * f1)
    F2 .= F + (Δℓ/2 * F_prime - RR0T * f2)
    M1 .= M - (Δℓ/2 * M_prime - RR0T * m1)
    M2 .= M + (Δℓ/2 * M_prime - RR0T * m2)

    # Nodal rotation angles
    θ[1] = rotation_angle(p1)
    θ[2] = rotation_angle(p2)

end

# Computes the modal states (generalized displacements, forces, strains, velocities and momenta) of the element at midpoint and at the nodes
function element_modal_states!(element,modalState,eigenvector,forceScaling)

    @unpack Δℓ,k,R0,R0_n1,R0_n2,R0T_n1,R0T_n2 = element.geometry
    @unpack C11,C12,C21,C22,I11,I12,I21,I22 = element.sectional
    @unpack DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω = element.assemblyIndices
    @unpack f1,f2,m1,m2 = element.stateCache.loadsResultants
    @unpack u,p,F,M,V,Ω,γ,κ,P,H = modalState.midpoint
    @unpack u1,u2,u1_b,u2_b,p1,p2,p1_b,p2_b,F1,F2,M1,M2,θ = modalState.nodal

    # Steady-state midpoint states, about which the modal states (perturbations) are defined
    γS = element.state.midpoint.γ
    κS = element.state.midpoint.κ
    RR0S = element.stateCache.rotationMatrices.RR0
    HTinvS = element.stateCache.rotationMatrices.HTinv
    R_p1S = element.stateCache.rotationMatrices.R_p1
    R_p2S = element.stateCache.rotationMatrices.R_p2
    R_p3S = element.stateCache.rotationMatrices.R_p3
    HTinv_p1S = element.stateCache.rotationMatrices.HTinv_p1
    HTinv_p2S = element.stateCache.rotationMatrices.HTinv_p2
    HTinv_p3S = element.stateCache.rotationMatrices.HTinv_p3

    # Midpoint states 
    u .= eigenvector[DOF_u]
    p .= eigenvector[DOF_p]
    F .= eigenvector[DOF_F] * forceScaling
    M .= eigenvector[DOF_M] * forceScaling
    V .= eigenvector[DOF_V]
    Ω .= eigenvector[DOF_Ω]

    # Midpoint complementary states
    γ .= C11 * F + C12 * M
    κ .= C21 * F + C22 * M
    P .= I11 * V + I12 * Ω
    H .= I21 * V + I22 * Ω

    # Rotation tensors
    R = rotation_tensor_WM(p)
    RR0 = R * R0
    RR0T = RR0'

    # Midpoint generalized displacements' derivatives w.r.t x1, linearized about the steady state
    u_prime = mul3(R_p1S,R_p2S,R_p3S,SVector(R0 * (a1 + γS))) * p + RR0S * γ
    p_prime = mul3(HTinv_p1S,HTinv_p2S,HTinv_p3S,SVector(R0 * κS)) * p + HTinvS * SVector(R0 * κ)

    # Midpoint deformed curvature vector, resolved in basis B
    K = k + κ

    # Midpoint generalized forces' derivatives w.r.t x1 (not linearized, approximate)
    F_prime = cross(Ω,P) - cross(K,F)
    M_prime = cross(Ω,H) + cross(V,P) - cross(a1+γ,F) - cross(K,M)

    # Nodal displacements and rotation parameters' vectors, resolved in basis A
    u1 .= u - Δℓ/2 * u_prime
    u2 .= u + Δℓ/2 * u_prime
    p1 .= p - Δℓ/2 * p_prime
    p2 .= p + Δℓ/2 * p_prime

    # Nodal rotation angles
    θ[1] = rotation_angle(p1)
    θ[2] = rotation_angle(p2)

    # Nodal displacements and rotation parameters' vectors, resolved in basis b
    u1_b .= R0T_n1 * u1
    u2_b .= R0T_n2 * u2
    p1_b .= R0T_n1 * p1
    p2_b .= R0T_n2 * p2

    # Nodal generalized internal forces, resolved in basis B
    F1 .= F - (Δℓ/2 * F_prime - RR0T * f1)
    F2 .= F + (Δℓ/2 * F_prime - RR0T * f2)
    M1 .= M - (Δℓ/2 * M_prime - RR0T * m1)
    M2 .= M + (Δℓ/2 * M_prime - RR0T * m2)

end