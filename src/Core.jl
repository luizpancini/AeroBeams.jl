"""
element_arrays!(problem::Problem,model::Model,element::Element)

Computes the elemental contributions to the system's arrays (residual, Jacobian, inertia)

# Arguments
- problem::Problem
- model::Model
- element::Element
"""
function element_arrays!(problem::Problem,model::Model,element::Element)

    ## Generalized velocities of basis b at element's midpoint, resolved in basis A
    # --------------------------------------------------------------------------
    element_velocities_basis_b!(model,element,problem.σ,problem.timeNow)

    ## States' rates
    # --------------------------------------------------------------------------
    # States' rates
    element_states_rates!(problem,element)
    
    ## Rotation variables
    # --------------------------------------------------------------------------
    element_rotation_variables!(problem,element)
    
    ## Distributed loads
    # --------------------------------------------------------------------------
    element_distributed_loads!(problem,model,element)
    
    ## Check intent
    # If the intent is getting the external forces array, we need to calculate the loads with actual states (because follower loads depend on it), and then reset element states to zero (because R(x) = J*x - F_ext(x), then F_ext = -R only if x = 0)
    if problem.getExternalForcesArray == true
        u,p,F,M,V,Ω,udot,pdot,Vdot,Ωdot = zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3)
        RR0 = element.R0
        @pack! element = u,p,F,M,V,Ω,udot,pdot,Vdot,Ωdot,RR0
    end
    
    ## Sectional strains, momenta and momenta's rates
    # --------------------------------------------------------------------------
    # Strains
    element_strains!(element)
    # Momenta
    element_momenta!(element)
    # Momenta rates
    if problem isa DynamicProblem
        element_momenta_rates!(element)
    end
    
    ## Residual array
    # --------------------------------------------------------------------------
    element_residual!(problem,model,element)
    if problem.getExternalForcesArray || problem.skipJacobianUpdate
        return
    end
    
    ## Follower distributed loads' derivatives (including aerodynamic loads) w.r.t. the rotation parameters (contributions to Jacobian matrix)
    # --------------------------------------------------------------------------
    distributed_loads_derivatives_rotation_parameters!(element)

    ## Aerodynamic derivatives (contributions to Jacobian and inertia matrices)
    # ------------------------------------------------------------------------- 
    aero_derivatives!(problem,model,element)

    ## Jacobian matrix
    # --------------------------------------------------------------------------
    element_jacobian!(problem,model,element)
     
end


"""
special_node_arrays!(problem::Problem,model::Model,specialNode::SpecialNode)

Computes the nodal contributions to the system's arrays (residual, Jacobian, inertia)

# Arguments
- problem::Problem
- model::Model
- specialNode::SpecialNode
"""
function special_node_arrays!(problem::Problem,model::Model,specialNode::SpecialNode)

    # If the intent is getting the external forces vector, we need to calculate
    # the loads with actual states (because follower loads depend on it), and
    # then reset nodal displacements to zero (because since R(x) = J*x - F_ext(x), then F_ext = -R only if x = 0)
    if problem.getExternalForcesArray
        specialNode.u,specialNode.p = zeros(3),zeros(3)
    end

    # Add spring loads, if applicable
    if specialNode.hasSprings
        spring_loads!(model,specialNode)
    end

    ## Residual
    # --------------------------------------------------------------------------
    special_node_residual!(problem,model,specialNode)
    if problem.getExternalForcesArray || problem.skipJacobianUpdate
        return
    end

    ## Follower loads' derivatives w.r.t. the rotation parameters (contributions to Jacobian matrix)
    # --------------------------------------------------------------------------
    special_node_follower_loads_derivatives_rotation_parameters!(problem,specialNode)

    ## Jacobian
    # --------------------------------------------------------------------------
    special_node_jacobian!(problem,model,specialNode)

end


"""
element_velocities_basis_b!(model::Model,element::Element,σ::Float64=1.0,timeNow::Number=0.0)

Computes the generalized velocities of basis b at the element's midpoint, resolved in basis A

# Arguments
- model::Model
- element::Element
- σ::Float64 = load factor
- timeNow::Number = current time
"""
function element_velocities_basis_b!(model::Model,element::Element,σ::Float64=1.0,timeNow::Number=0.0)

    @unpack R_AT,v_A,ω_A = model
    @unpack r = element

    # Translational velocity
    v = σ * R_AT * (v_A(timeNow) + cross(ω_A(timeNow),R_AT'*r))

    # Angular velocity
    ω = σ * R_AT * ω_A(timeNow)

    @pack! element = v,ω

    return v,ω

end


"""
element_accelerations_basis_b!(model::Model,element::Element,σ::Float64=1.0,timeNow::Number=0.0)

Computes the generalized accelerations of basis b at the element's midpoint, resolved in basis A

# Arguments
- model::Model
- element::Element
- σ::Float64 = load factor
- timeNow::Number = current time
"""
function element_accelerations_basis_b!(model::Model,element::Element,σ::Float64=1.0,timeNow::Number=0.0)

    @unpack R_AT,v_A,ω_A,vdot_A,ωdot_A = model
    @unpack r = element

    # Translational acceleration
    vdot = σ * R_AT * (vdot_A(timeNow) + cross(ωdot_A(timeNow),R_AT'*r) + cross(ω_A(timeNow),cross(ω_A(timeNow),R_AT'*r)) - cross(ω_A(timeNow),R_AT'*v_A(timeNow)))

    # Angular acceleration
    ωdot = σ * R_AT * ωdot_A(timeNow)

    @pack! element = vdot,ωdot

    return vdot,ωdot

end


"""
element_states!(problem::Problem,model::Model,element::Element)

Gets the states (generalized displacements, forces, velocities and aerodynamic) of the element

# Arguments
- problem::Problem
- model::Model
- element::Element
"""
function element_states!(problem::Problem,model::Model,element::Element)

    @unpack x = problem
    @unpack forceScaling = model
    @unpack states,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,DOF_χ = element
    @unpack rotationConstraint = element

    u = x[DOF_u]
    p = x[DOF_p]
    if !isnothing(rotationConstraint)
        @unpack DOF,masterElementGlobalID,value = rotationConstraint
        p[DOF] = model.elements[masterElementGlobalID].states.p[DOF] + value
    end
    F = x[DOF_F]*forceScaling
    M = x[DOF_M]*forceScaling
    V = x[DOF_V]
    Ω = x[DOF_Ω]
    χ = !isempty(DOF_χ) ? x[DOF_χ] : Vector{Float64}()

    @pack! element.states = u,p,F,M,V,Ω,χ

end


"""
element_states_rates!(problem::Problem,element::Element)

Computes the states' rates of the current element

# Arguments
- problem::Problem
- element::Element
"""
function element_states_rates!(problem::Problem,element::Element)

    # Skip for all but dynamic problems
    if !isa(problem,DynamicProblem)
        return 
    end

    # Unpack
    @unpack Δt = problem   
    @unpack udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv,χdotEquiv = element
    @unpack u,p,V,Ω,χ = element.states
    @unpack udot,pdot,Vdot,Ωdot,χdot = element.statesRates
    
    # Current rates
    if !isinf(Δt)
        udot = 2/Δt*u - udotEquiv
        pdot = 2/Δt*p - pdotEquiv
        Vdot = 2/Δt*V - VdotEquiv
        Ωdot = 2/Δt*Ω - ΩdotEquiv
        χdot = 2/Δt*χ - χdotEquiv
    end

    @pack! element.statesRates = udot,pdot,Vdot,Ωdot,χdot

end
 

"""
element_rotation_variables!(problem::Problem,element::Element)

Computes the rotation variables for the current element

# Arguments
- problem::Problem
- element::Element
"""
function element_rotation_variables!(problem::Problem,element::Element)

    @unpack R0 = element
    @unpack p = element.states
    @unpack pdot = element.statesRates

    ## Rotation tensors
    # --------------------------------------------------------------------------
    # Rotation tensor from basis b to basis B, resolved in basis A, and associated variables
    R,Θ,pNorm,λ,ps,ps1,ps2,ps3,ps0,υ,υ²,ps1s,ps2s,ps3s,ps1ps2,ps2ps3,ps1ps3 = rotation_tensor_WM(p)
    # Rotation tensor from basis A to basis B, resolved in basis A, and its transpose
    RR0 = R*R0 
    RR0T = Matrix(RR0')
    
    ## Function of tangent operator tensors
    # --------------------------------------------------------------------------
    HT = tangent_operator_transpose_WM(ps,ps0,υ²)
    HTinv = tangent_operator_transpose_inverse_WM(ps,ps0)
    
    ## Rotation tensor derivative w.r.t. extended rotation parameters
    # --------------------------------------------------------------------------
    R_p1,R_p2,R_p3,R_ps1,R_ps2,R_ps3,υ²_ps1,υ²_ps2,υ²_ps3,Θ_ps1,Θ_ps2,Θ_ps3,ps_p,ps1_p1,ps2_p1,ps3_p1,ps1_p2,ps2_p2,ps3_p2,ps1_p3,ps2_p3,ps3_p3 = rotation_tensor_derivatives_extended_parameters(p,pNorm,λ,ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²)
    
    ## Tangent operator tensors derivative w.r.t. extended rotation parameters
    # --------------------------------------------------------------------------
    HT_p1,HT_p2,HT_p3,HTinv_p1,HTinv_p2,HTinv_p3 = tangent_tensor_functions_derivatives_extended_parameters(HT,ps1,ps2,ps3,υ²,υ²_ps1,υ²_ps2,υ²_ps3,ps_p)
    
    ## Time derivatives for dynamic analyses
    # --------------------------------------------------------------------------
    if problem isa DynamicProblem
        # Time derivatives of rotation tensor from basis A to basis B, and scaled rotation parameters 
        Rdot,ps1dot,ps2dot,ps3dot = rotation_tensor_time_derivative(R_ps1,R_ps2,R_ps3,ps_p,pdot) 
        RdotR0 = Rdot*R0
                                                                              
        # Derivatives of time derivative of rotation tensor from basis A to basis B w.r.t extended rotation parameters
        Rdot_p1,Rdot_p2,Rdot_p3 = rotation_tensor_derivatives_time_extended_parameters(ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,Θ_ps1,Θ_ps2,Θ_ps3,υ,υ²,υ²_ps1,υ²_ps2,υ²_ps3,ps1dot,ps2dot,ps3dot,ps1_p1,ps2_p1,ps3_p1,ps1_p2,ps2_p2,ps3_p2,ps1_p3,ps2_p3,ps3_p3)
    else
        RdotR0,Rdot_p1,Rdot_p2,Rdot_p3 = zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3)
    end

    @pack! element = R,RR0,RR0T,RdotR0,HT,HTinv,R_p1,R_p2,R_p3,HT_p1,HT_p2,HT_p3,HTinv_p1,HTinv_p2,HTinv_p3,Rdot_p1,Rdot_p2,Rdot_p3
    
end


"""
element_distributed_loads!(problem::Problem,model::Model,element::Element)

Computes the nodal resultants from distributed loads on the current element, resolved in basis A

# Arguments
- problem::Problem
- model::Model
- element::Element
"""
function element_distributed_loads!(problem::Problem,model::Model,element::Element)

    @unpack σ = problem

    # Gravitational loads 
    f_g,m_g = gravitational_loads!(model,element,σ)

    # Externally applied distributed loads
    f1_d,f2_d,m1_d,m2_d = distributed_external_loads!(problem,element,σ)

    # Aerodynamic loads
    f1_χ,f2_χ,m1_χ,m2_χ = aerodynamic_loads!(problem,model,element)

    # Total nodal resultants from distributed loads
    f1 = f_g + f1_d + f1_χ
    f2 = f_g + f2_d + f2_χ
    m1 = m_g + m1_d + m1_χ
    m2 = m_g + m2_d + m2_χ

    @pack! element = f1,f2,m1,m2

end


"""
gravitational_loads!(model::Model,element::Element,σ::Float64)

Computes the nodal resultants from the distributed gravitational loads on the current element

# Arguments
- model::Model
- element::Element
- σ::Float64
"""
function gravitational_loads!(model::Model,element::Element,σ::Float64)

    @unpack gravityVector,R_AT = model
    @unpack Δℓ,μ,ηtilde,RR0,RR0T = element

    # Nodal distributed weight force and moment vectors, resolved in basis A
    if any(!iszero(gravityVector)) 
        f_g = σ * Δℓ/2 * μ * R_AT * gravityVector
        m_g = σ * RR0 * ηtilde * RR0T * f_g
    else
        f_g,m_g = zeros(3),zeros(3)
    end

    @pack! element = f_g,m_g

    return f_g,m_g
end


"""
distributed_external_loads!(problem::Problem,element::Element,σ::Float64)

Computes the nodal resultants from the externally applied distributed loads on the current element

# Arguments
- model::Model
- element::Element
- σ::Float64
"""
function distributed_external_loads!(problem::Problem,element::Element,σ::Float64)

    @unpack R,R0,RR0,f_A,m_A,f_b,m_b,ff_A,mf_A,ff_b,mf_b,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb = element

    # Initialize nodal resultants resolved in basis A
    f1_d,f2_d,m1_d,m2_d = zeros(3),zeros(3),zeros(3),zeros(3)

    # Initialize current values of nodal resultants from distributed follower loads
    ff1_A,ff2_A,mf1_A,mf2_A,ff1_b,ff2_b,mf1_b,mf2_b = zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3)

    # Add dead forces initially resolved in basis A
    if hasDistributedDeadForcesBasisA
        # Interpolate, if needed
        f1_A,f2_A = interpolate_distributed_loads(problem,σ*f_A)
        # Add
        f1_d += f1_A
        f2_d += f2_A
    end

    # Add dead moments initially resolved in basis A
    if hasDistributedDeadMomentsBasisA
        # Interpolate, if needed
        m1_A,m2_A = interpolate_distributed_loads(problem,σ*m_A)
        # Add
        m1_d += m1_A
        m2_d += m2_A
    end

    # Add dead forces initially resolved in basis b
    if hasDistributedDeadForcesBasisb
        # Interpolate, if needed
        f1_b,f2_b = interpolate_distributed_loads(problem,σ*f_b)
        # Add
        f1_d += R0 * f1_b
        f2_d += R0 * f2_b
    end

    # Add dead moments initially resolved in basis b
    if hasDistributedDeadMomentsBasisb
        # Interpolate, if needed
        m1_b,m2_b = interpolate_distributed_loads(problem,σ*m_b)
        # Add
        m1_d += R0 * m1_b
        m2_d += R0 * m2_b
    end

    # Add follower forces initially resolved in basis A
    if hasDistributedFollowerForcesBasisA
        # Interpolate, if needed
        ff1_A,ff2_A = interpolate_distributed_loads(problem,σ*ff_A)
        # Add
        f1_d += R * ff1_A
        f2_d += R * ff2_A
    end

    # Add follower moments initially resolved in basis A
    if hasDistributedFollowerMomentsBasisA
        # Interpolate, if needed
        mf1_A,mf2_A = interpolate_distributed_loads(problem,σ*mf_A)
        # Add
        m1_d += R * mf1_A
        m2_d += R * mf2_A
    end

    # Add follower forces initially resolved in basis b
    if hasDistributedFollowerForcesBasisb
        # Interpolate, if needed
        ff1_b,ff2_b = interpolate_distributed_loads(problem,σ*ff_b)
        # Add
        f1_d += RR0 * ff1_b
        f2_d += RR0 * ff2_b
    end

    # Add follower moments initially resolved in basis b
    if hasDistributedFollowerMomentsBasisb
        # Interpolate, if needed
        mf1_b,mf2_b = interpolate_distributed_loads(problem,σ*mf_b)
        # Add
        m1_d += RR0 * mf1_b
        m2_d += RR0 * mf2_b
    end

    @pack! element = ff1_A,ff2_A,mf1_A,mf2_A,ff1_b,ff2_b,mf1_b,mf2_b

    return f1_d,f2_d,m1_d,m2_d

end


"""
aerodynamic_loads!(problem::Problem,model::Model,element::Element)

Computes the nodal resultants from the aerodynamic loads on the current element

# Arguments
- problem::Problem
- model::Model
- element::Element
"""
function aerodynamic_loads!(problem::Problem,model::Model,element::Element)

    @unpack aero = element

    # Initialize 
    f1_χ,f2_χ,m1_χ,m2_χ = zeros(3),zeros(3),zeros(3),zeros(3)
    
    # Skip if there are no aerodynamic loads on the element
    if isnothing(aero)
        @pack! element = f1_χ,f2_χ,m1_χ,m2_χ
        return f1_χ,f2_χ,m1_χ,m2_χ
    end

    @unpack x,timeNow = problem
    @unpack DOF_δ = element
    @unpack V,Ω,χ = element.states
    @unpack Vdot,Ωdot = element.statesRates
    @unpack nTotalAeroStates,airfoil,δ,δMultiplier = aero
    δNow = !isempty(DOF_δ) ? x[DOF_δ]*δMultiplier : δ(timeNow)

    # Compute aerodynamic nodal resutants and state matrices
    f1_χ,f2_χ,m1_χ,m2_χ,A,B = aero_loads_core!(problem,model,element,V,Ω,χ,Vdot,Ωdot,δNow)

    @pack! element = f1_χ,f2_χ,m1_χ,m2_χ
    @pack! element.aero = A,B,δNow

    return f1_χ,f2_χ,m1_χ,m2_χ

end


"""
wrapper_aerodynamic_loads_from_states!(states,problem::Problem,model::Model,element::Element)

Wrapper function for aerodynamic loads using elemental states as inputs

# Arguments
- states
- problem::Problem
- model::Model
- element::Element
"""
function wrapper_aerodynamic_loads_from_states!(states,problem::Problem,model::Model,element::Element)

    # Unpack states and states' rates
    if isempty(element.DOF_δ)
        V,Ω,χ,δNow = states[1:3],states[4:6],states[7:end],element.aero.δ(problem.timeNow)
    else
        V,Ω,χ,δNow = states[1:3],states[4:6],states[7:end-1],states[end]
    end
    @unpack Vdot,Ωdot = element.statesRates

    # Nodal resutants and state matrices
    f1_χ,f2_χ,m1_χ,m2_χ,A,B = aero_loads_core!(problem,model,element,V,Ω,χ,Vdot,Ωdot,δNow)

    # Set outputs array
    out = vcat([f1_χ,f2_χ,m1_χ,m2_χ,vec(A),B]...)
    return out

end


"""
wrapper_aerodynamic_loads_from_states_rates!(statesRates,problem::Problem,model::Model,element::Element)

Wrapper function for aerodynamic loads using elemental states' rates as inputs

# Arguments
- statesRates
- problem::Problem
- model::Model
- element::Element
"""
function wrapper_aerodynamic_loads_from_states_rates!(statesRates,problem::Problem,model::Model,element::Element)

    # Unpack states and states' rates
    Vdot,Ωdot = statesRates[1:3],statesRates[4:6]
    @unpack V,Ω,χ = element.states
    δNow = !isempty(element.DOF_δ) ? problem.x[element.DOF_δ]*element.aero.δMultiplier : element.aero.δ(problem.timeNow)

    # Nodal resutants and state matrices
    f1_χ,f2_χ,m1_χ,m2_χ,A,B = aero_loads_core!(problem,model,element,V,Ω,χ,Vdot,Ωdot,δNow)

    # Set outputs array
    out = vcat([f1_χ,f2_χ,m1_χ,m2_χ,vec(A),B]...)
    return out

end


"""
aero_loads_core!(problem::Problem,model::Model,element::Element,V,Ω,χ,Vdot,Ωdot,δNow)

Computes the nodal resutants from aerodynamic loads and aerodynamic state matrices

# Arguments
- problem::Problem
- model::Model
- element::Element
- V
- Ω
- χ
- Vdot
- Ωdot
- δNow
"""
function aero_loads_core!(problem::Problem,model::Model,element::Element,V,Ω,χ,Vdot,Ωdot,δNow)

    # Steady aerodynamic kinematics
    aero_steady_kinematics!(element,V,Ω)

    # Skip if angle of attack is not defined
    if isnan(element.aero.flowAnglesAndRates.α)
        f1_χ,f2_χ,m1_χ,m2_χ,A,B = zeros(3),zeros(3),zeros(3),zeros(3),zeros(element.aero.nTotalAeroStates,element.aero.nTotalAeroStates),zeros(element.aero.nTotalAeroStates)
        return f1_χ,f2_χ,m1_χ,m2_χ,A,B
    end
    
    # Unsteady aerodynamic kinematics
    aero_unsteady_kinematics!(element,Vdot,Ωdot)

    # Non-dimensional flow parameters
    nondimensional_flow_parameters!(model,element)

    # Update airfoil parameters
    update_airfoil_parameters!(problem,element)

    # Local gust velocity
    # local_gust_velocity!(model,element,timeNow)

    # Flap deflection rates
    if !element.aero.δIsZero
        flap_deflection_rates!(problem,element)
    end

    # Effective (unsteady) angle of attack
    effective_angle_of_attack!(element,χ,δNow)

    # Aerodynamic coefficients 
    aero_coefficients!(problem,element,χ,δNow)

    # Aerodynamic state matrices 
    A,B = aero_state_matrices!(element,δNow)
    
    # Aerodynamic loads' nodal resultants
    f1_χ,f2_χ,m1_χ,m2_χ = aero_loads_resultants!(model,element)

    return f1_χ,f2_χ,m1_χ,m2_χ,A,B

end


"""
interpolate_distributed_loads(problem::Problem,loadArray::Array{Float64})

Interpolates the loads array at the current time

# Arguments
- problem::Problem
- loadArray::Array{Float64}
"""
function interpolate_distributed_loads(problem::Problem,loadArray::Array{Float64})

    # For steady problems, no interpolation is needed
    if !isa(problem,DynamicProblem)
        v1 = loadArray[1,:,1]
        v2 = loadArray[2,:,1]
        return v1, v2
    end

    # Unpack
    @unpack timeNow,timeBeginTimeStep,timeEndTimeStep,indexBeginTimeStep,indexEndTimeStep = problem

    # Check if the current time is between time steps
    if (timeBeginTimeStep < timeNow < timeEndTimeStep)
        # Initialize time vector and respective indices' range for interpolations
        timesBeginAndEnd = [timeBeginTimeStep,timeEndTimeStep]
        indicesBeginAndEnd = indexBeginTimeStep:indexEndTimeStep
        # Interpolate
        v1 = [interpolate(timesBeginAndEnd,loadArray[1,i,indicesBeginAndEnd],timeNow) for i in 1:3]
        v2 = [interpolate(timesBeginAndEnd,loadArray[2,i,indicesBeginAndEnd],timeNow) for i in 1:3]
    else
        # Set values at the begin of the time step
        v1 = loadArray[1,:,indexBeginTimeStep]
        v2 = loadArray[2,:,indexBeginTimeStep]
    end

    return v1, v2
end


"""
element_strains!(element::Element)

Computes the strains for the current element, resolved in basis B

# Arguments
- element::Element
"""
function element_strains!(element::Element)

    @unpack states,compStates,S = element
    @unpack F,M = states
    
    # Generalized sectional forces array
    forces = [F; M]
    
    # Sectional strains vector 
    strains = S*forces
    γ = strains[1:3] 
    κ = strains[4:6]

    @pack! element.compStates = γ,κ
    
end


"""
element_momenta!(element::Element)

Computes the momenta for the current element, resolved in basis B

# Arguments
- element::Element
"""
function element_momenta!(element::Element)

    @unpack states,compStates,I = element
    @unpack V,Ω = states

    # Generalized sectional velocities array
    velocities = [V; Ω]
    
    # Sectional linear and angular momenta vector 
    momenta = I*velocities
    P = momenta[1:3]
    H = momenta[4:6]

    @pack! element.compStates = P,H
    
end


"""
element_momenta_rates!(element::Element)

Computes the momenta rates for the current element, resolved in basis B

# Arguments
- element::Element
"""
function element_momenta_rates!(element::Element)

    @unpack statesRates,compStatesRates,I = element
    @unpack Vdot,Ωdot = statesRates

    # Generalized sectional accelerations array
    accelerations = [Vdot; Ωdot]
    
    # Sectional linear and angular momenta rates vector 
    momentaRates = I*accelerations
    Pdot = momentaRates[1:3]
    Hdot = momentaRates[4:6]

    @pack! element.compStatesRates = Pdot,Hdot
    
end


"""
element_residual!(problem::Problem,model::Model,element::Element)

Computes the contributions from the current element to the residual array

# Arguments
- problem::Problem
- model::Model
- element::Element
"""
function element_residual!(problem::Problem,model::Model,element::Element)

    @unpack residual = problem
    @unpack forceScaling = model
    @unpack Δℓ,R0,R0T,RR0,RR0T,RdotR0,HT,HTinv,v,ω,f1,f2,m1,m2,eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,eqs_FF1_sep,eqs_FF2_sep,eqs_FM1_sep,eqs_FM2_sep,eqs_Fχ,isSpecialNode1,isSpecialNode2,eqsNode1Set,eqsNode2Set,hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat = element
    @unpack u,p,F,M,V,Ω,χ = element.states
    @unpack γ,κ,P,H = element.compStates
    @unpack udot,pdot,χdot = element.statesRates
    @unpack Pdot,Hdot = element.compStatesRates

    ## Static terms
    # --------------------------------------------------------------------------
    # --- F_u --- #
    tmp = RR0 * F
    F_u1 = -tmp - f1
    F_u2 =  tmp - f2

    # --- F_p --- #
    tmp = RR0 * M
    tmp2 = Δℓ/2 * RR0 * cross(a1+γ,F)
    F_p1 = -tmp - m1 - tmp2
    F_p2 =  tmp - m2 - tmp2

    # --- F_F --- #
    tmp = Δℓ/2 * (RR0 * (a1+γ) - R0*a1)
    F_F1 =  u - tmp
    F_F2 = -u - tmp

    # --- F_M --- #
    tmp = Δℓ/2 * HTinv * R0 * κ
    F_M1 =  p - tmp
    F_M2 = -p - tmp

    ## Steady terms
    # ------------------------------------------------------------------------- 
    # --- F_u --- #
    tmp = Δℓ/2 * cross(ω,RR0*P)
    F_u1 += tmp
    F_u2 += tmp

    # --- F_p --- #
    tmp = Δℓ/2 * (cross(ω,RR0*H) + RR0 * cross(V,P))
    F_p1 += tmp
    F_p2 += tmp

    # --- F_V --- #
    F_V = RR0*V - v - cross(ω,u)
    
    # --- F_Ω --- #
    F_Ω = Ω - RR0T*ω

    # --- F_χ --- #
    if !isempty(eqs_Fχ)
        @unpack A,B = element.aero
        F_χ = -(A*χ+B)
    end

    ## Transient dynamic terms
    # -------------------------------------------------------------------------
    if problem isa DynamicProblem
        # --- F_u --- #
        tmp = Δℓ/2 * (RdotR0*P + RR0*Pdot)
        F_u1 += tmp
        F_u2 += tmp

        # --- F_p --- #    
        tmp = Δℓ/2 * (RdotR0*H + RR0*Hdot)
        F_p1 += tmp
        F_p2 += tmp

        # --- F_V --- #
        F_V += -udot
        
        # --- F_Ω --- #
        F_Ω += -R0T * HT * pdot

        # --- F_χ --- #
        if !isempty(eqs_Fχ)
            F_χ += χdot
        end
    end

    ## Insert element resultants from dynamic equilibrium and generalized displacements compatibility equations into the residual array
    # -------------------------------------------------------------------------

    # Set or add to residuals for the element's first node
    # --------------------------------------------------------------------------
    # The node's equations have not been set yet
    if !eqsNode1Set            
        # Set equilibrium equations' residual
        residual[eqs_Fu1] = F_u1/forceScaling
        residual[eqs_Fp1] = F_p1/forceScaling
        # Set compatibility equations' residual
        residual[eqs_FF1] = F_F1
        residual[eqs_FM1] = F_M1 
    # The node's equations have already been set    
    else                      
        # Add to existing equilibrium equations' residual (except for hinged degrees-of-freedom)
        residual[eqs_Fu1] += F_u1/forceScaling
        residual[eqs_Fp1] += notHingedNode1Mat*F_p1/forceScaling
        # Check if is a special node 
        # ----------------------------------------------------------------------
        # This is a standard node (shares compatibility equations)
        if !isSpecialNode1
            # Add to existing compatibility equations' residual / set moment equilibrium equations' residual for hinged degrees-of-freedom
            residual[eqs_FF1] += F_F1
            residual[eqs_FM1] += notHingedNode1Mat*F_M1 .+ hingedNode1Mat*F_p1/forceScaling 
        # This is a special node (has separate compatibility equations)    
        else
            # Set separate compatibility equations' residual  
            residual[eqs_FF1_sep] = F_F1
            residual[eqs_FM1_sep] = F_M1
        end
    end
    # Set or add to residuals for the element's last node
    # --------------------------------------------------------------------------
    # The node's equations have not been set yet
    if !eqsNode2Set          
        # Set equilibrium equations' residual
        residual[eqs_Fu2] = F_u2/forceScaling
        residual[eqs_Fp2] = F_p2/forceScaling
        # Set compatibility equations' residual (except for hinged degrees-of-freedom)
        residual[eqs_FF2] = F_F2
        residual[eqs_FM2] = notHingedNode2Mat*F_M2 
    # The node's equations have already been set    
    else                    
        # Add to existing equilibrium equations' residual
        residual[eqs_Fu2] += F_u2/forceScaling
        residual[eqs_Fp2] += F_p2/forceScaling
        # Check if is a special node 
        # ----------------------------------------------------------------------
        # This is a standard node (shares compatibility equations)
        if !isSpecialNode2     
            # Add to existing compatibility equations' residual 
            residual[eqs_FF2] += F_F2
            residual[eqs_FM2] += F_M2
        # This is a special node (has separate compatibility equations)    
        else              
            # Set separate compatibility equations' residual 
            residual[eqs_FF2_sep] = F_F2
            residual[eqs_FM2_sep] = F_M2
        end
    end

    # Generalized velocity-displacement residuals
    residual[eqs_FV] = F_V
    residual[eqs_FΩ] = F_Ω

    # Aerodynamic states residuals
    if !isempty(eqs_Fχ)
        residual[eqs_Fχ] = F_χ
    end

    @pack! problem = residual
end


"""
distributed_loads_derivatives_rotation_parameters!(element::Element)

Computes the derivatives of the distributed loads w.r.t the rotation parameters

# Arguments
- element::Element
"""
function distributed_loads_derivatives_rotation_parameters!(element::Element)

    # Derivatives of gravitational loads w.r.t. extended rotation parameters
    mg_p = gravitational_loads_derivatives_rotation_parameters(element)

    # Derivatives of externally applied follower distributed loads w.r.t. extended rotation parameters
    f1d_p,f2d_p,m1d_p,m2d_p = distributed_external_loads_derivatives_rotation_parameters(element)

    # Derivatives of aerodynamic loads w.r.t. extended rotation parameters
    f1χ_p,f2χ_p,m1χ_p,m2χ_p = aero_loads_derivatives_rotation_parameters(element)

    # Total nodal resultants' derivatives w.r.t. extended rotation parameters
    f1_p = f1χ_p + f1d_p
    m1_p = mg_p + m1χ_p + m1d_p
    f2_p = f2χ_p + f2d_p
    m2_p = mg_p + m2χ_p + m2d_p

    @pack! element = f1_p,f2_p,m1_p,m2_p

end


"""
gravitational_loads_derivatives_rotation_parameters(element::Element)

Computes the contributions of the gravitational loads to the Jacobian matrix

# Arguments
- element::Element
"""
function gravitational_loads_derivatives_rotation_parameters(element::Element)

    @unpack R,RR0,R0T,R_p1,R_p2,R_p3,f_g,m_g,ηtilde = element

    # Distributed weight moment vector derivative w.r.t. extended rotation parameters
    if any(!iszero(f_g))
        mg_p = mul3(R_p1,R_p2,R_p3,R'*m_g) + RR0*ηtilde*R0T*mul3(Matrix(R_p1'),Matrix(R_p2'),Matrix(R_p3'),f_g)
    else
        mg_p = zeros(3,3)
    end

    return mg_p

end


"""
distributed_external_loads_derivatives_rotation_parameters(element::Element)

Computes the contributions of the externally applied distributed loads to the Jacobian matrix

# Arguments
- element::Element
"""
function distributed_external_loads_derivatives_rotation_parameters(element::Element)

    @unpack ff1_A,ff2_A,mf1_A,mf2_A,ff1_b,ff2_b,mf1_b,mf2_b,R_p1,R_p2,R_p3,R0,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb = element

    # Initialize derivatives of externally applied distributed loads w.r.t. extended rotation parameters
    f1d_p,f2d_p,m1d_p,m2d_p = zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3)

    # Contributions from follower forces initially resolved in basis A
    if hasDistributedFollowerForcesBasisA
        f1d_p += mul3(R_p1,R_p2,R_p3,ff1_A)
        f2d_p += mul3(R_p1,R_p2,R_p3,ff2_A)
    end

    # Contributions from follower moments initially resolved in basis A
    if hasDistributedFollowerMomentsBasisA
        m1d_p += mul3(R_p1,R_p2,R_p3,mf1_A)
        m2d_p += mul3(R_p1,R_p2,R_p3,mf2_A)
    end

    # Contributions from follower forces initially resolved in basis b
    if hasDistributedFollowerForcesBasisb
        f1d_p += mul3(R_p1,R_p2,R_p3,R0*ff1_b)
        f2d_p += mul3(R_p1,R_p2,R_p3,R0*ff2_b)
    end

    # Contributions from follower moements initially resolved in basis b
    if hasDistributedFollowerMomentsBasisb
        m1d_p += mul3(R_p1,R_p2,R_p3,R0*mf1_b)
        m2d_p += mul3(R_p1,R_p2,R_p3,R0*mf2_b)
    end

    return f1d_p,f2d_p,m1d_p,m2d_p

end


"""
aero_loads_derivatives_rotation_parameters(element::Element)

Computes the derivatives of the aerodynamic loads w.r.t the extended rotation parameters

# Arguments
- element::Element
"""
function aero_loads_derivatives_rotation_parameters(element::Element)

    @unpack aero = element

    # Skip if there are no aero loads
    if isnothing(aero)
        f1χ_p,f2χ_p,m1χ_p,m2χ_p = zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3)
        return f1χ_p,f2χ_p,m1χ_p,m2χ_p
    end

    @unpack R_p1,R_p2,R_p3 = element
    @unpack RwR0,F = aero

    f1χ_p = mul3(R_p1,R_p2,R_p3,RwR0*F[1:3])
    m1χ_p = mul3(R_p1,R_p2,R_p3,RwR0*F[4:6])
    f2χ_p = mul3(R_p1,R_p2,R_p3,RwR0*F[7:9])
    m2χ_p = mul3(R_p1,R_p2,R_p3,RwR0*F[10:12])

    return f1χ_p,f2χ_p,m1χ_p,m2χ_p

end


"""
aero_derivatives!(problem::Problem,model::Model,element::Element)

Computes the derivatives of the aerodynamic loads and aerodynamic states w.r.t remaining elemental states and their rates (excludes rotation parameters, already computed)

# Arguments
- problem::Problem
- model::Model
- element::Element
"""
function aero_derivatives!(problem::Problem,model::Model,element::Element)

    @unpack aero = element

    # Skip if there are no aero loads, or angle of attack is undefined, or the convergence rate in a dynamic problem is higher than the required minimum
    if isnothing(aero) || isnan(aero.flowAnglesAndRates.α) || (problem isa DynamicProblem && problem.systemSolver.convRate > problem.systemSolver.minConvRateAeroJacUpdate)
        return
    end

    # Unpack element data
    @unpack DOF_δ = element
    @unpack V,Ω,χ = element.states
    @unpack nTotalAeroStates,derivationMethod,A,δMultiplier = aero

    # Set input states for aero loads wrapper functions
    states = isempty(DOF_δ) ? vcat([V,Ω,χ]...) : vcat([V,Ω,χ,problem.x[DOF_δ]*δMultiplier]...)

    # Save complementary variables of dynamic stall model
    if typeof(element.aero.solver) in [BLi]
        BLcompVars = deepcopy(element.aero.BLcompVars)
    end

    # Get derivatives of aerodynamic loads and state matrices w.r.t. states
    if typeof(derivationMethod) == AD
        derivatives = ForwardDiff.jacobian(x -> wrapper_aerodynamic_loads_from_states!(x,problem,model,element), states)
    elseif typeof(derivationMethod) == FD
        derivatives = FiniteDifferences.jacobian(derivationMethod.method, x -> wrapper_aerodynamic_loads_from_states!(x,problem,model,element), states)
        derivatives = derivatives[1]
    end

    # Deal with possible NaN values
    NaNind = LinearIndices(derivatives)[findall(isnan,derivatives)]
    while !isempty(NaNind)
        if NaNind[1] == 1
            derivatives[1] = 0.0
            popat!(NaNind,1)
        end
        derivatives[NaNind] .= derivatives[NaNind .- 1]
        NaNind = LinearIndices(derivatives)[findall(isnan,derivatives)]
    end

    # Extract derivatives
    f1χ_V = derivatives[1:3,1:3]
    f2χ_V = derivatives[4:6,1:3]
    m1χ_V = derivatives[7:9,1:3]
    m2χ_V = derivatives[10:12,1:3]
    f1χ_Ω = derivatives[1:3,4:6]
    f2χ_Ω = derivatives[4:6,4:6]
    m1χ_Ω = derivatives[7:9,4:6]
    m2χ_Ω = derivatives[10:12,4:6]
    f1χ_χ = derivatives[1:3,7:6+nTotalAeroStates]
    f2χ_χ = derivatives[4:6,7:6+nTotalAeroStates]
    m1χ_χ = derivatives[7:9,7:6+nTotalAeroStates]
    m2χ_χ = derivatives[10:12,7:6+nTotalAeroStates]
    A_V = reshape(derivatives[13:12+nTotalAeroStates^2,1:3],(nTotalAeroStates,nTotalAeroStates,3))
    A_Ω = reshape(derivatives[13:12+nTotalAeroStates^2,4:6],(nTotalAeroStates,nTotalAeroStates,3))
    B_V = derivatives[13+nTotalAeroStates^2:end,1:3]
    B_Ω = derivatives[13+nTotalAeroStates^2:end,4:6]
    if !isempty(DOF_δ)
        f1χ_δ = derivatives[1:3,end]
        f2χ_δ = derivatives[4:6,end]
        m1χ_δ = derivatives[7:9,end]
        m2χ_δ = derivatives[10:12,end]
    else
        f1χ_δ,f2χ_δ,m1χ_δ,m2χ_δ = zeros(3),zeros(3),zeros(3),zeros(3)
    end

    # Set derivatives of aerodynamic states w.r.t. states
    F_χ_χ = -A
    F_χ_V,F_χ_Ω = zeros(nTotalAeroStates,3),zeros(nTotalAeroStates,3)
    if !isempty(χ)
        for i=1:3
            F_χ_V[:,i] = -(A_V[:,:,i]*χ+B_V[:,i])
            F_χ_Ω[:,i] = -(A_Ω[:,:,i]*χ+B_Ω[:,i])
        end
    end

    # Reset complementary variables of dynamic stall model
    if typeof(element.aero.solver) in [BLi]
        element.aero.BLcompVars = BLcompVars
    end

    @pack! element.aero = f1χ_V,f2χ_V,f1χ_Ω,f2χ_Ω,m1χ_V,m2χ_V,m1χ_Ω,m2χ_Ω,f1χ_χ,f2χ_χ,m1χ_χ,m2χ_χ,f1χ_δ,f2χ_δ,m1χ_δ,m2χ_δ,F_χ_V,F_χ_Ω,F_χ_χ

    # Skip if not an EigenProblem
    if !isa(problem,EigenProblem)
        return
    end

    # Unpack element data
    @unpack Vdot,Ωdot = element.statesRates

    # Set input states' rates for aero loads wrapper function
    statesRates = vcat([Vdot,Ωdot]...)

    # Get derivatives of aerodynamic loads and state matrices w.r.t. states' rates
    if typeof(derivationMethod) == AD
        derivatives = ForwardDiff.jacobian(x -> wrapper_aerodynamic_loads_from_states_rates!(x,problem,model,element), statesRates)
    elseif typeof(derivationMethod) == FD
        derivatives = FiniteDifferences.jacobian(derivationMethod.method, x -> wrapper_aerodynamic_loads_from_states_rates!(x,problem,model,element), statesRates)
        derivatives = derivatives[1]
    end

    # Deal with possible NaN values
    NaNind = LinearIndices(derivatives)[findall(isnan,derivatives)]
    while !isempty(NaNind)
        if NaNind[1] == 1
            derivatives[1] = 0.0
            popat!(NaNind,1)
        end
        derivatives[NaNind] .= derivatives[NaNind .- 1]
        NaNind = LinearIndices(derivatives)[findall(isnan,derivatives)]
    end

    # Extract derivatives
    f1χ_Vdot = derivatives[1:3,1:3]
    f2χ_Vdot = derivatives[4:6,1:3]
    m1χ_Vdot = derivatives[7:9,1:3]
    m2χ_Vdot = derivatives[10:12,1:3]
    f1χ_Ωdot = derivatives[1:3,4:6]
    f2χ_Ωdot = derivatives[4:6,4:6]
    m1χ_Ωdot = derivatives[7:9,4:6]
    m2χ_Ωdot = derivatives[10:12,4:6]
    A_Vdot = reshape(derivatives[13:12+nTotalAeroStates^2,1:3],(nTotalAeroStates,nTotalAeroStates,3))
    A_Ωdot = reshape(derivatives[13:12+nTotalAeroStates^2,4:6],(nTotalAeroStates,nTotalAeroStates,3))
    B_Vdot = derivatives[13+nTotalAeroStates^2:end,1:3]
    B_Ωdot = derivatives[13+nTotalAeroStates^2:end,4:6]

    # Set derivatives of aerodynamic states w.r.t. states' rates
    F_χ_Vdot,F_χ_Ωdot = zeros(nTotalAeroStates,3),zeros(nTotalAeroStates,3)
    F_χ_χdot = Matrix(1.0*LinearAlgebra.I,nTotalAeroStates,nTotalAeroStates)
    if !isempty(χ)
        for i=1:3
            F_χ_Vdot[:,i] = -(A_Vdot[:,:,i]*χ+B_Vdot[:,i])
            F_χ_Ωdot[:,i] = -(A_Ωdot[:,:,i]*χ+B_Ωdot[:,i])
        end
    end

    # Reset complementary variables of dynamic stall model
    if typeof(element.aero.solver) in [BLi]
        element.aero.BLcompVars = BLcompVars
    end

    @pack! element.aero = f1χ_Vdot,f2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Vdot,m2χ_Vdot,m1χ_Ωdot,m2χ_Ωdot,F_χ_Vdot,F_χ_Ωdot,F_χ_χdot

end


"""
element_jacobian!(problem::Problem,model::Model,element::Element)

Computes the contributions from the current element to the Jacobian matrix

# Arguments
- problem::Problem
- model::Model
- element::Element
"""
function element_jacobian!(problem::Problem,model::Model,element::Element)

    @unpack jacobian = problem
    @unpack forceScaling = model
    @unpack Δℓ,S_11,S_12,S_21,S_22,I_11,I_12,I_21,I_22,ω,R0,R0T,RR0,RdotR0,HT,HTinv,R_p1,R_p2,R_p3,HT_p1,HT_p2,HT_p3,HTinv_p1,HTinv_p2,HTinv_p3,Rdot_p1,Rdot_p2,Rdot_p3,f1_p,f2_p,m1_p,m2_p,eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FF1_sep,eqs_FF2_sep,eqs_FM1_sep,eqs_FM2_sep,eqs_FV,eqs_FΩ,eqs_Fχ,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,DOF_χ,DOF_δ,eqsNode1Set,eqsNode2Set,isSpecialNode1,isSpecialNode2,hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat,aero = element
    @unpack F,M,V,Ω = element.states
    @unpack pdot,Vdot,Ωdot = element.statesRates
    @unpack γ,κ,P,H = element.compStates
    @unpack Pdot,Hdot = element.compStatesRates

    ## Static terms
    # --------------------------------------------------------------------------
    # Repeatedly used variables
    F_tilde = tilde(F)
        
    # --- F_u --- #
    # F_u_p
    tmp = mul3(R_p1,R_p2,R_p3,R0*F)
    F_u1_p = -tmp - f1_p
    F_u2_p =  tmp - f2_p
    # F_u_F
    F_u1_F = -RR0
    F_u2_F =  RR0

    # --- F_p --- #
    # F_p_p
    tmp1 = mul3(R_p1,R_p2,R_p3,R0*M)
    tmp2 = Δℓ/2 * mul3(R_p1,R_p2,R_p3,R0*cross(a1+γ,F))
    F_p1_p = -tmp1 - m1_p - tmp2
    F_p2_p =  tmp1 - m2_p - tmp2
    # F_p_F
    F_p1_F = -Δℓ/2 * RR0 * (tilde(a1+γ)-F_tilde*S_11)
    F_p2_F = F_p1_F
    # F_p_M
    tmp = Δℓ/2 * RR0 * F_tilde * S_12
    F_p1_M = tmp - RR0
    F_p2_M = tmp + RR0

    # --- F_F --- #
    # F_F_u
    F_F1_u = I3
    F_F2_u = -F_F1_u
    # F_F_p
    F_F1_p = -mul3(R_p1,R_p2,R_p3,Δℓ/2*R0*(a1+γ))
    F_F2_p = F_F1_p
    # F_F_F
    F_F1_F = -Δℓ/2 * RR0 * S_11
    F_F2_F = F_F1_F
    # F_F_M
    F_F1_M = -Δℓ/2 * RR0 * S_12
    F_F2_M = F_F1_M

    # --- F_M --- #
    # F_M_p
    tmp = mul3(HTinv_p1,HTinv_p2,HTinv_p3,Δℓ/2*R0*κ)
    F_M1_p =  I3 - tmp
    F_M2_p = -I3 - tmp
    # F_M_F
    tmp = -Δℓ/2 * HTinv * R0
    F_M1_F = tmp * S_21
    F_M2_F = F_M1_F
    # F_M_M
    F_M1_M = tmp * S_22
    F_M2_M = F_M1_M

    ## Steady-state terms
    # --------------------------------------------------------------------------
    # Repeatedly used variables
    ω_tilde = tilde(ω)
    ω_tilde_RR0 = ω_tilde*RR0
    V_tilde = tilde(V)
    R0P = R0*P
    R0H = R0*H
    
    # --- F_u --- #
    # F_u_p
    tmp = Δℓ/2 * ω_tilde * mul3(R_p1,R_p2,R_p3,R0P)
    F_u1_p += tmp
    F_u2_p += tmp 
    # F_u_V
    tmp = Δℓ/2 * ω_tilde_RR0 * I_11
    F_u1_V = tmp
    F_u2_V = tmp 
    # F_u_Ω
    tmp = Δℓ/2 * ω_tilde_RR0 * I_12
    F_u1_Ω = tmp
    F_u2_Ω = tmp 
    
    # --- F_p --- # 
    # F_p_p
    tmp = Δℓ/2 * (ω_tilde*mul3(R_p1,R_p2,R_p3,R0H) + mul3(R_p1,R_p2,R_p3,R0*cross(V,P)))
    F_p1_p += tmp
    F_p2_p += tmp  
    # F_p_V
    tmp = Δℓ/2 * (ω_tilde_RR0*I_21 + RR0*(V_tilde*I_11-tilde(P)))
    F_p1_V = tmp 
    F_p2_V = tmp    
    # F_p_Ω
    tmp = Δℓ/2 * (ω_tilde_RR0*I_22 + RR0*V_tilde*I_12)
    F_p1_Ω = tmp 
    F_p2_Ω = tmp 
    
    # --- F_V --- #    
    # F_V_u
    F_V_u = -ω_tilde   
    # F_V_p
    F_V_p = mul3(R_p1,R_p2,R_p3,R0*V) 
    # F_V_V
    F_V_V = RR0
    
    # --- F_Ω --- #   
    # F_Ω_p
    F_Ω_p = -R0T * mul3(Matrix(R_p1'),Matrix(R_p2'),Matrix(R_p3'),ω)  
    # F_Ω_Ω
    F_Ω_Ω = I3  

    ## Transient dynamic terms
    # --------------------------------------------------------------------------
    if problem isa DynamicProblem

        @unpack Δt = problem
        
        # Repeatedly used variables
        RdotR0_plus_2oΔtRR0 = RdotR0 + 2/Δt*RR0
        
        # --- F_u --- #
        # F_u_p
        tmp = Δℓ/2 * (mul3(Rdot_p1,Rdot_p2,Rdot_p3,R0P) + mul3(R_p1,R_p2,R_p3,R0*Pdot))
        F_u1_p += tmp
        F_u2_p += tmp
        # F_u_V
        tmp = Δℓ/2 * RdotR0_plus_2oΔtRR0 * I_11
        F_u1_V += tmp
        F_u2_V += tmp
        # F_u_Ω
        tmp = Δℓ/2 * RdotR0_plus_2oΔtRR0 * I_12
        F_u1_Ω += tmp
        F_u2_Ω += tmp  

        # --- F_p --- #
        # F_p_p
        tmp = Δℓ/2 * (mul3(Rdot_p1,Rdot_p2,Rdot_p3,R0H) + mul3(R_p1,R_p2,R_p3,R0*Hdot))
        F_p1_p += tmp
        F_p2_p += tmp
        # F_p_V
        tmp = Δℓ/2 * RdotR0_plus_2oΔtRR0 * I_21
        F_p1_V += tmp
        F_p2_V += tmp
        # F_p_Ω
        tmp = Δℓ/2 * RdotR0_plus_2oΔtRR0 * I_22
        F_p1_Ω += tmp
        F_p2_Ω += tmp
        
        # --- F_V --- #
        # F_V_u
        F_V_u += -2/Δt * I3

        # --- F_Ω --- #
        # F_Ω_p
        F_Ω_p += -R0T * (mul3(HT_p1,HT_p2,HT_p3,pdot) + 2/Δt*HT)

    end

    ## Aerodynamic terms 
    # --------------------------------------------------------------------------
    if !isnothing(aero)
        @unpack nTotalAeroStates,f1χ_V,f2χ_V,f1χ_Ω,f2χ_Ω,m1χ_V,m2χ_V,m1χ_Ω,m2χ_Ω,f1χ_χ,f2χ_χ,m1χ_χ,m2χ_χ,F_χ_V,F_χ_Ω,F_χ_χ = element.aero
        # F_u_V
        F_u1_V -= f1χ_V
        F_u2_V -= f2χ_V 
        # F_u_Ω
        F_u1_Ω -= f1χ_Ω
        F_u2_Ω -= f2χ_Ω
        # F_p_V
        F_p1_V -= m1χ_V
        F_p2_V -= m2χ_V   
        # F_p_Ω
        F_p1_Ω -= m1χ_Ω
        F_p2_Ω -= m2χ_Ω
        # F_χ_χ 
        if problem isa DynamicProblem && nTotalAeroStates > 0
            F_χ_χ += Matrix(2/Δt*LinearAlgebra.I,nTotalAeroStates,nTotalAeroStates)
        end
    end


    ## Insert element resultants into the Jacobian matrix
    # --------------------------------------------------------------------------

    ## Dynamic equilibrium and generalized displacements compatibility equations' contributions
    # --------------------------------------------------------------------------
    # ---------- Node 1 ----------- #
    # Equilibrium equations' Jacobian entries for the element's first node
    # --- F_u1 --- #
    jacobian[eqs_Fu1,DOF_p] = F_u1_p/forceScaling
    jacobian[eqs_Fu1,DOF_F] = F_u1_F
    jacobian[eqs_Fu1,DOF_V] = F_u1_V/forceScaling
    jacobian[eqs_Fu1,DOF_Ω] = F_u1_Ω/forceScaling
    # --- F_p1 --- #
    jacobian[eqs_Fp1,DOF_p] = notHingedNode1Mat*F_p1_p/forceScaling
    jacobian[eqs_Fp1,DOF_F] = notHingedNode1Mat*F_p1_F
    jacobian[eqs_Fp1,DOF_M] = notHingedNode1Mat*F_p1_M
    jacobian[eqs_Fp1,DOF_V] = F_p1_V/forceScaling
    jacobian[eqs_Fp1,DOF_Ω] = F_p1_Ω/forceScaling
    # If the first compatibility equations' Jacobians have already been set for this special node, use node's separate compatibility equations
    if eqsNode1Set && isSpecialNode1 
        eqs_FF1 = eqs_FF1_sep                         
        eqs_FM1 = eqs_FM1_sep
    end
    # Compatibility equations' Jacobian entries for the element's first node
    # --- F_F1 --- #
    jacobian[eqs_FF1,DOF_u] = F_F1_u
    jacobian[eqs_FF1,DOF_p] = F_F1_p
    jacobian[eqs_FF1,DOF_F] = F_F1_F*forceScaling
    jacobian[eqs_FF1,DOF_M] = F_F1_M*forceScaling
    # --- F_M1 --- #
    jacobian[eqs_FM1,DOF_p] = notHingedNode1Mat*F_M1_p .+ hingedNode1Mat*F_p1_p/forceScaling
    jacobian[eqs_FM1,DOF_F] = notHingedNode1Mat*F_M1_F*forceScaling .+ hingedNode1Mat*F_p1_F
    jacobian[eqs_FM1,DOF_M] = notHingedNode1Mat*F_M1_M*forceScaling  .+ hingedNode1Mat*F_p1_M

    # ---------- Node 2 ----------- #
    # Equilibrium equations' Jacobian entries for the element's second node
    # --- F_u2 --- #
    jacobian[eqs_Fu2,DOF_p] = F_u2_p/forceScaling
    jacobian[eqs_Fu2,DOF_F] = F_u2_F
    jacobian[eqs_Fu2,DOF_V] = F_u2_V/forceScaling
    jacobian[eqs_Fu2,DOF_Ω] = F_u2_Ω/forceScaling
    # --- F_p2 --- #
    jacobian[eqs_Fp2,DOF_p] = F_p2_p/forceScaling
    jacobian[eqs_Fp2,DOF_F] = F_p2_F
    jacobian[eqs_Fp2,DOF_M] = F_p2_M
    jacobian[eqs_Fp2,DOF_V] = F_p2_V/forceScaling
    jacobian[eqs_Fp2,DOF_Ω] = F_p2_Ω/forceScaling
    # If the first compatibility equations' Jacobians have already been set for this special node, use node's separate compatibility equations
    if eqsNode2Set && isSpecialNode2  
        eqs_FF2 = eqs_FF2_sep                         
        eqs_FM2 = eqs_FM2_sep
    end
    # Compatibility equations' Jacobian entries for the element's second node
    # --- F_F2 --- #
    jacobian[eqs_FF2,DOF_u] = F_F2_u
    jacobian[eqs_FF2,DOF_p] = F_F2_p
    jacobian[eqs_FF2,DOF_F] = F_F2_F*forceScaling
    jacobian[eqs_FF2,DOF_M] = F_F2_M*forceScaling
    # --- F_M2 --- #
    jacobian[eqs_FM2,DOF_p] = notHingedNode2Mat*F_M2_p
    jacobian[eqs_FM2,DOF_F] = notHingedNode2Mat*F_M2_F*forceScaling
    jacobian[eqs_FM2,DOF_M] = notHingedNode2Mat*F_M2_M*forceScaling

    ## Generalized velocities' Jacobian contributions
    # --------------------------------------------------------------------------
    # --- F_V --- #
    jacobian[eqs_FV,DOF_u] = F_V_u
    jacobian[eqs_FV,DOF_p] = F_V_p
    jacobian[eqs_FV,DOF_V] = F_V_V
    # --- F_Ω --- #
    jacobian[eqs_FΩ,DOF_p] = F_Ω_p
    jacobian[eqs_FΩ,DOF_Ω] = F_Ω_Ω

    ## Aerodynamic states' Jacobians
    # --------------------------------------------------------------------------
    if !isempty(DOF_χ)
        # --- F_u --- #
        jacobian[eqs_Fu1,DOF_χ] = -f1χ_χ/forceScaling
        jacobian[eqs_Fu2,DOF_χ] = -f2χ_χ/forceScaling
        # --- F_p --- #
        jacobian[eqs_Fp1,DOF_χ] = -m1χ_χ/forceScaling
        jacobian[eqs_Fp2,DOF_χ] = -m2χ_χ/forceScaling
        # --- F_χ --- #
        jacobian[eqs_Fχ,DOF_V] = F_χ_V
        jacobian[eqs_Fχ,DOF_Ω] = F_χ_Ω
        jacobian[eqs_Fχ,DOF_χ] = F_χ_χ
    end

    ## Flap trim state's Jacobians
    if !isempty(DOF_δ)
        @unpack f1χ_δ,f2χ_δ,m1χ_δ,m2χ_δ = element.aero
        # --- F_u --- #
        jacobian[eqs_Fu1,DOF_δ] = -f1χ_δ/forceScaling
        jacobian[eqs_Fu2,DOF_δ] = -f2χ_δ/forceScaling
        # --- F_p --- #
        jacobian[eqs_Fp1,DOF_δ] = -m1χ_δ/forceScaling
        jacobian[eqs_Fp2,DOF_δ] = -m2χ_δ/forceScaling
    end

    @pack! problem = jacobian

end


"""
element_inertia!(problem::Problem,model::Model,element::Element)

Computes the contributions from the current element to the inertia matrix

# Arguments
- problem::Problem
- model::Model
- element::Element
"""
function element_inertia!(problem::Problem,model::Model,element::Element)

    @unpack inertia = problem
    @unpack forceScaling = model
    @unpack Δℓ,I_11,I_12,I_21,I_22,R0,R0T,RR0,HT,HTinv,ω,R_p1,R_p2,R_p3,eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FV,eqs_FΩ,eqs_Fχ,DOF_u,DOF_p,DOF_V,DOF_Ω,DOF_χ,eqsNode1Set,eqsNode2Set,isSpecialNode1,isSpecialNode2,aero = element
    @unpack V,Ω = element.states
    @unpack Vdot,Ωdot = element.statesRates
    @unpack P,H = element.compStates

    ## Structural terms
    # --------------------------------------------------------------------------
    # --- F_u --- #
    # F_u_pdot
    F_u1_pdot = Δℓ/2 * mul3(R_p1,R_p2,R_p3,R0*P)
    F_u2_pdot = F_u1_pdot
    # F_u_Vdot
    tmp = Δℓ/2 * RR0 * I_11
    F_u1_Vdot = tmp
    F_u2_Vdot = tmp
    # F_u_Ωdot
    tmp = Δℓ/2 * RR0 * I_12
    F_u1_Ωdot = tmp
    F_u2_Ωdot = tmp

    # --- F_p --- #
    # F_p_pdot
    F_p1_pdot = Δℓ/2 * mul3(R_p1,R_p2,R_p3,R0*H)
    F_p2_pdot = F_p1_pdot
    # F_p_Vdot
    tmp = Δℓ/2 * RR0 * I_21
    F_p1_Vdot = tmp
    F_p2_Vdot = tmp
    # F_p_Ωdot
    tmp = Δℓ/2 * RR0 * I_22
    F_p1_Ωdot = tmp
    F_p2_Ωdot = tmp

    # --- F_V --- #
    # F_V_udot
    F_V_udot = -I3

    # --- F_Ω --- #
    # F_Ω_pdot
    F_Ω_pdot = -R0T * HT

    ## Aerodynamic loads terms
    # --------------------------------------------------------------------------
    if !isnothing(aero)
        @unpack f1χ_Vdot,f2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Vdot,m2χ_Vdot,m1χ_Ωdot,m2χ_Ωdot = element.aero
        # F_u_Vdot
        F_u1_Vdot -= f1χ_Vdot
        F_u2_Vdot -= f2χ_Vdot
        # F_u_Ωdot
        F_u1_Ωdot -= f1χ_Ωdot
        F_u2_Ωdot -= f2χ_Ωdot
        # F_p_Vdot
        F_p1_Vdot -= m1χ_Vdot
        F_p2_Vdot -= m2χ_Vdot
        # F_p_Ωdot
        F_p1_Ωdot -= m1χ_Ωdot
        F_p2_Ωdot -= m2χ_Ωdot
    end

    ## Insert element resultants into the inertia matrix
    # --------------------------------------------------------------------------
    # First node
    inertia[eqs_Fu1,DOF_p] = F_u1_pdot/forceScaling
    inertia[eqs_Fu1,DOF_V] = F_u1_Vdot/forceScaling
    inertia[eqs_Fu1,DOF_Ω] = F_u1_Ωdot/forceScaling
    inertia[eqs_Fp1,DOF_p] = F_p1_pdot/forceScaling
    inertia[eqs_Fp1,DOF_V] = F_p1_Vdot/forceScaling
    inertia[eqs_Fp1,DOF_Ω] = F_p1_Ωdot/forceScaling

    # Second node
    inertia[eqs_Fu2,DOF_p] = F_u2_pdot/forceScaling
    inertia[eqs_Fu2,DOF_V] = F_u2_Vdot/forceScaling
    inertia[eqs_Fu2,DOF_Ω] = F_u2_Ωdot/forceScaling
    inertia[eqs_Fp2,DOF_p] = F_p2_pdot/forceScaling
    inertia[eqs_Fp2,DOF_V] = F_p2_Vdot/forceScaling
    inertia[eqs_Fp2,DOF_Ω] = F_p2_Ωdot/forceScaling

    # Element's midpoint
    inertia[eqs_FV,DOF_u] = F_V_udot
    inertia[eqs_FΩ,DOF_p] = F_Ω_pdot
    if !isempty(eqs_Fχ)
        @unpack F_χ_Vdot,F_χ_Ωdot,F_χ_χdot = element.aero
        inertia[eqs_Fχ,DOF_V] = F_χ_Vdot
        inertia[eqs_Fχ,DOF_Ω] = F_χ_Ωdot
        inertia[eqs_Fχ,DOF_χ] = F_χ_χdot
    end

    @pack! problem = inertia

end


"""
element_nodal_states!(element::Element)

Updates the nodal states of the element

# Arguments
- element::Element
"""
function element_nodal_states!(element::Element)

    @unpack Δℓ,k,R0,RR0,RR0T,HTinv,f1,f2,m1,m2,R0_n1,R0_n2,R0T_n1,R0T_n2,nodalStates = element
    @unpack u,p,F,M,V,Ω = element.states
    @unpack γ,κ,P,H = element.compStates
    @unpack Pdot,Hdot = element.compStatesRates

    # Midpoint generalized displacements' derivatives with respect to x1
    u_prime = RR0*(a1+γ)-R0*a1
    p_prime = HTinv*R0*κ

    # Midpoint deformed curvature vector, resolved in basis B
    K = k+κ

    # Midpoint generalized forces' derivatives with respect to x1
    F_prime = Pdot+cross(Ω,P)-cross(K,F)
    M_prime = Hdot+cross(Ω,H)+cross(V,P)-cross(a1+γ,F)-cross(K,M)

    # Nodal displacements and rotation parameters' vectors, resolved in basis A
    u_n1 = u - Δℓ/2 * u_prime
    u_n2 = u + Δℓ/2 * u_prime
    p_n1 = p - Δℓ/2 * p_prime
    p_n2 = p + Δℓ/2 * p_prime

    # Nodal rotation angles
    θ_n1 = rotation_angle(p_n1)
    θ_n2 = rotation_angle(p_n2)

    # Nodal displacements and rotation parameters' vectors, resolved in basis b
    u_n1_b = R0T_n1*u_n1
    u_n2_b = R0T_n2*u_n2
    p_n1_b = R0T_n1*p_n1
    p_n2_b = R0T_n2*p_n2

    # Nodal generalized internal forces, resolved in basis B
    F_n1 = F - (Δℓ/2*F_prime - RR0T*f1)
    F_n2 = F + (Δℓ/2*F_prime - RR0T*f2)
    M_n1 = M - (Δℓ/2*M_prime - RR0T*m1)
    M_n2 = M + (Δℓ/2*M_prime - RR0T*m2)
 
    @pack! element.nodalStates = u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2

end


"""
special_node_states!(problem::Problem,model::Model,specialNode::SpecialNode)

Computes the nodal states

# Arguments
- problem::Problem
- model::Model
- specialNode::SpecialNode
"""
function special_node_states!(problem::Problem,model::Model,specialNode::SpecialNode)

    @unpack x,σ = problem
    @unpack uIsPrescribed,pIsPrescribed,trimIsPrescribed,BCs,DOF_uF,DOF_pM,DOF_trimLoads = specialNode

    # Check if node is not BC'ed - all nodal states are generalized displacements, externally applied forces are null
    if isempty(BCs)
        u = x[DOF_uF]
        p = x[DOF_pM]
        F = zeros(3)
        M = zeros(3)
        @pack! specialNode = u,p,F,M
        return
    end

    @unpack forceScaling = model

    # Initialize states  
    u,p,F,M = zeros(3),zeros(3),zeros(3),zeros(3)

    # Loop boundary conditions
    for BC in BCs

        @unpack currentValue,isLoad,isFollower,isTrim,R0_n = BC

        ## Generalized displacements
        # ----------------------------------------------------------------------
        # Loop directions and set displacements and rotations, resolved in basis A, as either nodal states or prescribed values
        for i=1:3
            if uIsPrescribed[i] && !isLoad[i]
                # u[i] is a prescribed value
                u[i] = σ*currentValue[i]
            elseif !uIsPrescribed[i] && isLoad[i]
                # u[i] is a nodal state      
                u[i] = x[DOF_uF[i]]                 
            end
            if pIsPrescribed[i] && !isLoad[i+3]
                # p[i] is a prescribed value
                p[i] = σ*currentValue[i+3]        
            elseif !pIsPrescribed[i] && isLoad[i+3]
                # p[i] is a nodal state 
                p[i] = x[DOF_pM[i]]       
            end
        end

        ## Prescribed values of generalized forces, resolved in basis A
        # ----------------------------------------------------------------------
        # Nodal rotation tensor from basis b to basis B 
        R,_ = rotation_tensor_WM(p) 

        # Initialize prescribed nodal values (hatted values) for the loads of current BC
        F̂,M̂ = zeros(3),zeros(3)

        # Initialize trim nodal values
        Ftrim,Mtrim = zeros(3),zeros(3)
        
        # Loop orthogonal directions
        for i=1:3
            # Add non-trim forces
            if isLoad[i] && !isTrim[i]
                # Add dead forces (dead values are constant when resolved in basis A)
                F̂ += σ*currentValue[i]*!isFollower[i]*I3[:,i]
                # Add follower forces (follower values are constant when resolved in basis B - transform to basis A)
                F̂ += σ*currentValue[i]*isFollower[i]*R[:,i]
            end
            # Add trim forces
            if isTrim[i]
                # Add dead trim forces
                Ftrim += x[DOF_trimLoads[i]]*!isFollower[i]*I3[:,i]*forceScaling
                # Add follower trim forces
                Ftrim += x[DOF_trimLoads[i]]*isFollower[i]*R[:,i]*forceScaling
            end
            # Add non-trim moments
            if isLoad[i+3] && !isTrim[i+3]
                # Add dead moments (dead values are constant when resolved in basis A)
                M̂ += σ*currentValue[i+3]*!isFollower[i+3]*I3[:,i]
                # Add follower moments (follower values are constant when resolved in basis B - transform to basis A)
                M̂ += σ*currentValue[i+3]*isFollower[i+3]*R[:,i]
            end
            # Add trim moments
            if isTrim[i+3]
                # Add dead trim moments
                Mtrim += x[DOF_trimLoads[i+3]]*!isFollower[i+3]*I3[:,i]*forceScaling
                # Add follower trim moments
                Mtrim += x[DOF_trimLoads[i+3]]*isFollower[i+3]*R[:,i]*forceScaling
            end
        end
    
        ## Generalized forces ------------------------------------------------------------------------
        # Loop directions and set forces and moments, resolved in basis A, as either nodal states or prescribed values
        for i=1:3
            if isLoad[i] && !isTrim[i]
                # F[i] is a prescribed value
                F[i] += F̂[i]  
            elseif isTrim[i]
                # F[i] is a trim value
                F[i] = Ftrim[i]                   
            elseif !trimIsPrescribed[i]
                # F[i] is a nodal state
                F[i] = x[DOF_uF[i]]*forceScaling
            end
            if isLoad[i+3] && !isTrim[i+3]
                # M[i] is a prescribed value
                M[i] += M̂[i] 
            elseif isTrim[i+3] 
                # M[i] is a trim value
                M[i] = Mtrim[i]                   
            elseif !trimIsPrescribed[i+3]
                # M[i] is a nodal state
                M[i] = x[DOF_pM[i]]*forceScaling   
            end
        end
    end

    @pack! specialNode = u,p,F,M
end


"""
spring_loads!(model::Model,specialNode::SpecialNode)

Adds the contributions of attached springs to the special node's loads

# Arguments
- model::Model
- specialNode::SpecialNode
"""
function spring_loads!(model::Model,specialNode::SpecialNode)

    @unpack u,p,F,M,springs,globalID = specialNode

    # Loop springs
    for spring in springs
        @unpack hasDoubleAttachment,Ku,Kp = spring
        # Double attachment
        if hasDoubleAttachment
            @unpack nodesSpecialIDs,nodesGlobalIDs = spring
            # Node of other attachment
            otherNode = globalID == nodesGlobalIDs[1] ? model.specialNodes[nodesSpecialIDs[2]] : model.specialNodes[nodesSpecialIDs[1]]
            # Displacement and rotation of that node
            uOtherNode = otherNode.u
            pOtherNode = otherNode.p
            # Add spring loads, resolved in basis A
            F -= Ku * (u - uOtherNode)
            M -= Kp * (p - pOtherNode)
        # Single attachment    
        else
            # Add spring loads, resolved in basis A
            F -= Ku * u
            M -= Kp * p
        end
    end

    @pack! specialNode = F,M
end


"""
special_node_residual!(problem::Problem,model::Model,specialNode::SpecialNode)

Computes the contributions from the current node to the residual array

# Arguments
- problem::Problem
- model::Model
- specialNode::SpecialNode
"""
function special_node_residual!(problem::Problem,model::Model,specialNode::SpecialNode)

    @unpack residual = problem
    @unpack forceScaling = model
    @unpack ζonElements,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,eqs_FF_sep,eqs_FM_sep,u,p,F,M = specialNode

    # Loop connected elements of the node
    for e in eachindex(ζonElements)
        # Check if the node has separate compatibility equations
        if isempty(eqs_FF_sep[e])   
            # Add to node's equilibrium and compatibility equations' residuals
            residual[eqs_Fu[e]] -= F/forceScaling
            residual[eqs_Fp[e]] -= M/forceScaling
            residual[eqs_FF[e]] += ζonElements[e]*u
            residual[eqs_FM[e]] += ζonElements[e]*p
        else                    
            # Add to node's separate compatibility equations
            residual[eqs_FF_sep[e]] += ζonElements[e]*u
            residual[eqs_FM_sep[e]] += ζonElements[e]*p
        end
    end

    @pack! problem = residual

end


"""
special_node_follower_loads_derivatives_rotation_parameters!(problem::Problem,specialNode::SpecialNode)

Computes the contributions from the nodal follower loads to the Jacobian matrix

# Arguments
- problem::Problem
- specialNode::SpecialNode
"""
function special_node_follower_loads_derivatives_rotation_parameters!(problem::Problem,specialNode::SpecialNode)

    @unpack x,σ = problem
    @unpack forceScaling = problem.model
    @unpack BCs,p,DOF_trimLoads = specialNode

    # Initialize
    F_p,M_p = zeros(3,3),zeros(3,3)

    # Check if the node is BC'ed or if there are any follower loads
    followers = vcat([BC.isFollower for BC in BCs]...)
    if isempty(BCs) || !any(followers)
        @pack! specialNode = F_p,M_p
        return 
    end

    # Get nodal rotation tensor's derivatives w.r.t. the extended rotation parameters
    R_p1,R_p2,R_p3 = rotation_tensor_derivatives_extended_parameters(p)

    # Loop boundary conditions
    for BC in BCs

        @unpack currentValue,isLoad,isFollower,isTrim,R0_n = BC

        # Skip if there are no follower loads
        if !any(isFollower)
            continue
        end

        # Nodal follower load components in the undeformed configuration, resolved in basis A
        F₀,M₀ = zeros(3),zeros(3)
        for i=1:3
            # Forces
            if isTrim[i] && isFollower[i]  
                # Trim variable 
                F₀[i] = x[DOF_trimLoads[i]]*forceScaling
            elseif isFollower[i]           
                # Not a trim variable
                F₀[i] = σ*currentValue[i]
            end
            # Moments
            if isTrim[i+3] && isFollower[i+3]
                # Trim variable 
                M₀[i] = x[DOF_trimLoads[i+3]]*forceScaling
            elseif isFollower[i+3]       
                # Not a trim variable
                M₀[i] = σ*currentValue[i]
            end
        end

        # Add to derivatives of the follower loads w.r.t. the extended rotation parameters 
        F_p += mul3(R_p1,R_p2,R_p3,F₀)
        M_p += mul3(R_p1,R_p2,R_p3,M₀)
    end

    @pack! specialNode = F_p,M_p

end


"""
special_node_jacobian!(problem::Problem,model::Model,specialNode::SpecialNode)

Computes the contributions from the current node to the Jacobian matrix

# Arguments
- problem::Problem
- model::Model
- specialNode::SpecialNode
"""
function special_node_jacobian!(problem::Problem,model::Model,specialNode::SpecialNode)

    @unpack jacobian = problem
    @unpack forceScaling = model
    @unpack globalID,ζonElements,BCs,uIsPrescribed,pIsPrescribed,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,eqs_FF_sep,eqs_FM_sep,DOF_uF,DOF_pM,DOF_trimLoads,u,p,F,M,F_p,M_p,hasSprings,springs = specialNode

    # Check if the node is BC'ed
    if !isempty(BCs)
        # Loop applied BCs
        for BC in BCs
            @unpack isLoad,isTrim = BC
            jacobian = update_special_node_jacobian!(jacobian,forceScaling,F_p,M_p,ζonElements,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,eqs_FF_sep,eqs_FM_sep,DOF_uF,DOF_pM,DOF_trimLoads,uIsPrescribed,pIsPrescribed,isLoad,isTrim)
        end
    else
        jacobian = update_special_node_jacobian!(jacobian,forceScaling,F_p,M_p,ζonElements,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,eqs_FF_sep,eqs_FM_sep,DOF_uF,DOF_pM,DOF_trimLoads,uIsPrescribed,pIsPrescribed)
    end

    # Add spring loads' contributions, if applicable
    if hasSprings
        jacobian = spring_loads_jacobians!(model,jacobian,forceScaling,globalID,eqs_Fu,eqs_Fp,DOF_uF,DOF_pM,springs)
    end

    @pack! problem = jacobian

end


"""
update_special_node_jacobian!(jacobian,forceScaling,F_p,M_p,ζonElements,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,eqs_FF_sep,eqs_FM_sep,DOF_uF,DOF_pM,DOF_trimLoads,uIsPrescribed,pIsPrescribed,isLoad=trues(6),isTrim=falses(6))

Updates the Jacobian matrix with the contributions from the current BC (if any) on the current node 

# Arguments
- problem::Problem
- model::Model
- specialNode::SpecialNode
"""
function update_special_node_jacobian!(jacobian,forceScaling,F_p,M_p,ζonElements,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,eqs_FF_sep,eqs_FM_sep,DOF_uF,DOF_pM,DOF_trimLoads,uIsPrescribed,pIsPrescribed,isLoad=trues(6),isTrim=falses(6))

    # Loop connected elements of the node
    for e in eachindex(ζonElements)
        # If the node does not have separate compatibility equations
        if isempty(eqs_FF_sep[e])
            # Loop directions and set Jacobian components for the node's equilibrium and compatibility equations
            for i=1:3
                # Forces/displacements DOFs
                if isLoad[i] && !uIsPrescribed[i]
                    # If F[i] is prescribed, u[i] is the nodal state: Set compatibility equation's Jacobian component [d(-/+u)/d(u) = -/+1 = ζonElements[e]] 
                    jacobian[eqs_FF[e][i],DOF_uF[i]] = ζonElements[e] 
                    for j=1:3 
                        # Loop over node's moments/rotations DOFs
                        if isLoad[j+3] && !pIsPrescribed[j]
                            # If M(j) is prescribed, p(j) is the nodal state: Set equilibrium equation's Jacobian component (follower load component)
                            jacobian[eqs_Fu[e][i],DOF_pM[j]] = -F_p[i,j]/forceScaling  
                        end
                    end
                else                                                        
                    # If u[i] is prescribed, F[i] is the nodal state: Set equilibrium equation's Jacobian component [d(-F)/d(F) = -1]
                    jacobian[eqs_Fu[e][i],DOF_uF[i]] = -1   
                end
                if isTrim[i]
                    # If F[i] is also a trim variable: Set equilibrium equation's Jacobian component [d(-F)/d(F) = -1]
                    jacobian[eqs_Fu[e][i],DOF_trimLoads[i]] = -1
                end 
                # Moments/rotations DOFs
                if isLoad[i+3] && !pIsPrescribed[i]
                    # If M[i] is prescribed, p[i] is the nodal state: Set compatibility equation's Jacobian component [d(-/+p)/d(p) = -/+1 = ζonElements[e]] 
                    jacobian[eqs_FM[e][i],DOF_pM[i]] = ζonElements[e]
                    for j=1:3
                        # Loop over node's moments/rotations DOFs
                        if isLoad[j+3] && !pIsPrescribed[j]
                            # If M(j) is prescribed, p(j) is the nodal state: Add equilibrium equation's Jacobian component (follower load component)
                            jacobian[eqs_Fp[e][i],DOF_pM[j]] = -M_p[i,j]/forceScaling   
                        end
                    end
                else                                                        
                    # If p[i] is prescribed, M[i] is the nodal state: Set equilibrium equation's Jacobian component [d(-M)/d(M) = -1]
                    jacobian[eqs_Fp[e][i],DOF_pM[i]] = -1
                end
                if isTrim[i+3]
                    # If M[i] is a trim variable: Set equilibrium equation's Jacobian component [d(-M)/d(M) = -1]
                    jacobian[eqs_Fp[e][i],DOF_trimLoads[i+3]] = -1
                end
            end  
        # Else, set Jacobian components for the node's separate compatibility equations   
        else                                                                        # Loop over node's forces/displacements DOFs
            for i=1:3
                # If F[i] is prescribed, u[i] is the nodal state: Set separate compatibility equation's Jacobian component [d(-/+u)/d(u) = -/+1 = ζonElements[e]]              
                if isLoad[i] && !uIsPrescribed[i]                         
                    jacobian[eqs_FF_sep[e][i],DOF_uF[i]] = ζonElements[e]  
                end
            end
            # Loop over node's moments/rotations DOFs
            for i=1:3
                # If M[i] is prescribed, p[i] is the nodal state: Set separate compatibility equation's Jacobian component [d(-/+p)/d(p) = -/+1 = ζonElements[e]]
                if isLoad[i+3] && !pIsPrescribed[i]
                    jacobian[eqs_FM_sep[e][i],DOF_pM[i]] = ζonElements[e]
                end
            end
        end
    end

    return jacobian
end


"""
spring_loads_jacobians!(model,jacobian,forceScaling,globalID,eqs_Fu,eqs_Fp,DOF_uF,DOF_pM,springs)

Adds the contributions of the spring loads to the Jacobian matrix

# Arguments
- model
- jacobian
- forceScaling
- globalID
- eqs_Fu
- eqs_Fp
- DOF_uF
- DOF_pM
- springs
"""
function spring_loads_jacobians!(model,jacobian,forceScaling,globalID,eqs_Fu,eqs_Fp,DOF_uF,DOF_pM,springs)

    # Loop springs
    for spring in springs
        @unpack hasDoubleAttachment,Ku,Kp = spring
        # Add translational spring Jacobian: d(-F_spring)/d(u) = Ku/forceScaling
        jacobian[eqs_Fu[1],DOF_uF] += Ku/forceScaling
        # Add rotational spring Jacobian: d(-M_spring)/d(p) = Kp/forceScaling
        jacobian[eqs_Fp[1],DOF_pM] += Kp/forceScaling
        # If doubly-attached, add contributions from other node's DOFs
        if hasDoubleAttachment
            @unpack nodesSpecialIDs,nodesGlobalIDs = spring
            # Node of other attachment
            otherNode = globalID == nodesGlobalIDs[1] ? model.specialNodes[nodesSpecialIDs[2]] : model.specialNodes[nodesSpecialIDs[1]]
            # Add translational spring Jacobian: d(-F_spring)/d(uOtherNode) = -Ku/forceScaling
            jacobian[eqs_Fu[1],otherNode.DOF_uF] -= Ku/forceScaling
            # Add rotational spring Jacobian: d(-M_spring)/d(pOtherNode) = -Kp/forceScaling
            jacobian[eqs_Fp[1],otherNode.DOF_pM] -= Kp/forceScaling
        end
    end

    return jacobian
end


"""
element_modal_states(element::Element,eigenvector::Vector{T},forceScaling::Float64)

Computes the modal states (generalized displacements, forces, strains, velocities and momenta) of the element at midpoint and at the nodes

# Arguments
- element::Element
- eigenvector::Vector{T}
- forceScaling::Float64
"""
function element_modal_states(element::Element,eigenvector::Vector{T},forceScaling::Float64) where T<:Union{Float64,ComplexF64}

    @unpack DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,Δℓ,k,R0,S,I,f1,f2,m1,m2,R0_n1,R0_n2,R0T_n1,R0T_n2 = element

    # Midpoint states 
    u = eigenvector[DOF_u]
    p = eigenvector[DOF_p]
    F = eigenvector[DOF_F]*forceScaling
    M = eigenvector[DOF_M]*forceScaling
    V = eigenvector[DOF_V]
    Ω = eigenvector[DOF_Ω]

    # Midpoint complementary states
    strains = S*[F; M]
    momenta = I*[V; Ω]
    γ = strains[1:3] 
    κ = strains[4:6]
    P = momenta[1:3]
    H = momenta[4:6]

    # Rotation tensors
    R,_ = rotation_tensor_WM(p) 
    RR0 = R*R0
    RR0T = Matrix(RR0')
    HTinv = tangent_operator_transpose_inverse_WM(p)

    # Midpoint generalized displacements' derivatives with respect to x1
    u_prime = RR0*(a1+γ)-R0*a1
    p_prime = HTinv*R0*κ

    # Midpoint deformed curvature vector, resolved in basis B
    K = k+κ

    # Midpoint generalized forces' derivatives with respect to x1
    F_prime = cross(Ω,P)-cross(K,F)
    M_prime = cross(Ω,H)+cross(V,P)-cross(a1+γ,F)-cross(K,M)

    # Nodal displacements and rotation parameters' vectors, resolved in basis A
    u_n1 = u - Δℓ/2 * u_prime
    u_n2 = u + Δℓ/2 * u_prime
    p_n1 = p - Δℓ/2 * p_prime
    p_n2 = p + Δℓ/2 * p_prime

    # Nodal rotation angles
    θ_n1 = rotation_angle(p_n1)
    θ_n2 = rotation_angle(p_n2)

    # Nodal displacements and rotation parameters' vectors, resolved in basis b
    u_n1_b = R0T_n1*u_n1
    u_n2_b = R0T_n2*u_n2
    p_n1_b = R0T_n1*p_n1
    p_n2_b = R0T_n2*p_n2

    # Nodal generalized internal forces, resolved in basis B
    F_n1 = F - (Δℓ/2*F_prime - RR0T*f1)
    F_n2 = F + (Δℓ/2*F_prime - RR0T*f2)
    M_n1 = M - (Δℓ/2*M_prime - RR0T*m1)
    M_n2 = M + (Δℓ/2*M_prime - RR0T*m2)

    return u,p,F,M,V,Ω,γ,κ,P,H,u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2
end


"""
update_states!(problem::Problem)

Updates the elemental and nodal states

# Arguments
- problem::Problem
"""
function update_states!(problem::Problem)

    @unpack model = problem
    @unpack elements = model

    for element in elements

        ## Generalized velocities of basis b at element's midpoint, resolved in basis A
        # ----------------------------------------------------------------------
        # Velocities
        element_velocities_basis_b!(model,element,problem.σ,problem.timeNow)
        # Accelerations
        element_accelerations_basis_b!(model,element,problem.σ,problem.timeNow)

        ## States and states' rates
        # ----------------------------------------------------------------------
        # States
        element_states!(problem,model,element)
        # States' rates
        element_states_rates!(problem,element)
        
        ## Rotation variables
        # ----------------------------------------------------------------------
        element_rotation_variables!(problem,element)
        
        ## Distributed loads
        # ----------------------------------------------------------------------
        element_distributed_loads!(problem,model,element)
        
        ## Sectional strains, momenta and momenta's rates
        # ----------------------------------------------------------------------
        # Strains
        element_strains!(element)
        # Momenta
        element_momenta!(element)
        # Momenta rates
        element_momenta_rates!(element)

        ## Update nodal states of the element (mere outputs)
        # ----------------------------------------------------------------------
        element_nodal_states!(element)
    end

end