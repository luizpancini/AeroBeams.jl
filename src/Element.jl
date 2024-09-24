#
# mutable struct ElementalStates

#     ElementalStates composite type

# Fields
# - `u::Vector{T}` = displacements resolved in basis A
# - `p::Vector{T}` = rotation parameters resolved in basis A
# - `F::Vector{T}` = sectional forces resolved in basis B
# - `M::Vector{T}` = sectional moments resolved in basis B
# - `V::Vector{T}` = sectional velocities resolved in basis B
# - `Ω::Vector{T}` = sectional angular velocities resolved in basis B
# - `χ::Vector{T}` = aerodynamic states
#
mutable struct ElementalStates{T<:Union{Float64,ComplexF64}}
    
    # Fields
    u::Vector{T}
    p::Vector{T}
    F::Vector{T}
    M::Vector{T}
    V::Vector{T}
    Ω::Vector{T}
    χ::Vector{T}

    # Constructor
    function ElementalStates{T}() where T<:Union{Float64,ComplexF64}

        u,p,F,M,V,Ω,χ = zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(0)

        return new{T}(u,p,F,M,V,Ω,χ)
    end
end


#
# mutable struct ComplementaryElementalStates

#     ComplementaryElementalStates composite type

# Fields
# - `γ::Vector{T}` = force strains resolved in basis b
# - `κ::Vector{T}` = moment strains resolved in basis b
# - `P::Vector{T}` = linear momenta resolved in basis B
# - `H::Vector{T}` = angular momenta resolved in basis B
#
mutable struct ComplementaryElementalStates{T<:Union{Float64,ComplexF64}}
    
    # Fields
    γ::Vector{T}
    κ::Vector{T}
    P::Vector{T}
    H::Vector{T}

    # Constructor
    function ComplementaryElementalStates{T}() where T<:Union{Float64,ComplexF64}

        γ,κ,P,H = zeros(3),zeros(3),zeros(3),zeros(3)

        return new{T}(γ,κ,P,H)
    end
end
    

#
# mutable struct ElementalStatesRates

#     ElementalStatesRates composite type

# Fields
# - `udot::Vector{Float64}` = displacements time derivative
# - `pdot::Vector{Float64}` = rotation parameters time derivative
# - `Vdot::Vector{Float64}` = sectional velocities time derivative
# - `Ωdot::Vector{Float64}` = sectional angular velocities time derivative
# - `χdot::Vector{Float64}` = aerodynamic states time derivative
#
mutable struct ElementalStatesRates
    
    # Fields
    udot::Vector{Float64}
    pdot::Vector{Float64}
    Vdot::Vector{Float64}
    Ωdot::Vector{Float64}
    χdot::Vector{Float64}

    # Constructor
    function ElementalStatesRates() 

        udot,pdot,Vdot,Ωdot = zeros(3),zeros(3),zeros(3),zeros(3)

        χdot = zeros(0)

        return new(udot,pdot,Vdot,Ωdot,χdot)
    end
end


#
# mutable struct ComplementaryElementalStatesRates

#     ComplementaryElementalStatesRates composite type

# Fields
# - `Pdot::Vector{Float64}` = linear momenta time derivative
# - `Hdot::Vector{Float64}` = angular momenta time derivative
#
mutable struct ComplementaryElementalStatesRates
    
    # Fields
    Pdot::Vector{Float64}
    Hdot::Vector{Float64}

    # Constructor
    function ComplementaryElementalStatesRates() 

        Pdot,Hdot = zeros(3),zeros(3)

        return new(Pdot,Hdot)
    end
end


#
# mutable struct NodalStates

#     NodalStates composite type

# Fields
# - `u_n1::Vector{T}` = displacements of node 1 resolved in basis A
# - `u_n2::Vector{T}` = displacements of node 2 resolved in basis A
# - `p_n1::Vector{T}` = rotation parameters of node 1 resolved in basis A
# - `p_n2::Vector{T}` = rotation parameters of node 2 resolved in basis A
# - `u_n1_b::Vector{T}` = displacements of node 1 resolved in basis b
# - `u_n2_b::Vector{T}` = displacements of node 2 resolved in basis b
# - `p_n1_b::Vector{T}` = rotation parameters of node 1 resolved in basis b
# - `p_n2_b::Vector{T}` = rotation parameters of node 2 resolved in basis b
# - `F_n1::Vector{T}` = sectional forces of node 1 resolved in basis B
# - `F_n2::Vector{T}` = sectional forces of node 2 resolved in basis B
# - `M_n1::Vector{T}` = sectional moments of node 1 resolved in basis B
# - `M_n2::Vector{T}` = sectional moments of node 2 resolved in basis B
# - `θ_n1::T` = rotation angle of node 1
# - `θ_n2::T` = rotation angle of node 2
#
mutable struct NodalStates{T<:Union{Float64,ComplexF64}}
    
    # Fields
    u_n1::Vector{T}
    u_n2::Vector{T}
    p_n1::Vector{T}
    p_n2::Vector{T}
    u_n1_b::Vector{T}
    u_n2_b::Vector{T}
    p_n1_b::Vector{T}
    p_n2_b::Vector{T}
    F_n1::Vector{T}
    F_n2::Vector{T}
    M_n1::Vector{T}
    M_n2::Vector{T}
    θ_n1::T
    θ_n2::T

    # Constructor
    function NodalStates{T}() where T<:Union{Float64,ComplexF64}

        u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2 = zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),0.0,0.0

        return new{T}(u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2)
    end
end


#
# mutable struct Resultants

#     Resultants composite type

# Fields
# - `F_u1::Vector{Float64}` = first node translational equilibrium equation resultant
# - `F_u2::Vector{Float64}` = second node translational equilibrium equation resultant
# - `F_p1::Vector{Float64}` = first node rotational equilibrium equation resultant
# - `F_p2::Vector{Float64}` = second node rotational equilibrium equation resultant
# - `F_F1::Vector{Float64}` = first node translational compatibility equation resultant
# - `F_F2::Vector{Float64}` = second node translational compatibility equation resultant
# - `F_M1::Vector{Float64}` = first node rotational compatibility equation resultant
# - `F_M2::Vector{Float64}` = second node rotational compatibility equation resultant
# - `F_V::Vector{Float64}` = midpoint linear velocity-displacement compatibility equation resultant
# - `F_Ω::Vector{Float64}` = midpoint angular velocity-displacement compatibility equation resultant
mutable struct Resultants
    
    # Fields
    F_u1::Vector{Float64}
    F_u2::Vector{Float64}
    F_p1::Vector{Float64}
    F_p2::Vector{Float64}
    F_F1::Vector{Float64}
    F_F2::Vector{Float64}
    F_M1::Vector{Float64}
    F_M2::Vector{Float64}
    F_V::Vector{Float64}
    F_Ω::Vector{Float64}
    F_χ::Vector{Float64}

    # Constructor
    function Resultants(nTotalAeroStates)

        F_u1 = zeros(3)
        F_u2 = zeros(3)
        F_p1 = zeros(3)
        F_p2 = zeros(3)
        F_F1 = zeros(3)
        F_F2 = zeros(3)
        F_M1 = zeros(3)
        F_M2 = zeros(3)
        F_V = zeros(3)
        F_Ω = zeros(3)
        F_χ = zeros(nTotalAeroStates)

        return new(F_u1,F_u2,F_p1,F_p2,F_F1,F_F2,F_M1,F_M2,F_V,F_Ω,F_χ)
    end
end


#
# mutable struct Jacobians

#     Jacobians composite type

# Notes
# These are the partial derivatives (Jacobians) of the finite element equations w.r.t. the states and states' rates (only the structural ones)
mutable struct Jacobians
    
    # Fields
    F_u1_p::Vector{Float64}
    F_u2_p::Vector{Float64}
    F_u1_F::Vector{Float64}
    F_u2_F::Vector{Float64}
    F_u1_V::Vector{Float64}
    F_u2_V::Vector{Float64}
    F_u1_Ω::Vector{Float64}
    F_u2_Ω::Vector{Float64}
    F_p1_p::Vector{Float64}
    F_p2_p::Vector{Float64}
    F_p1_F::Vector{Float64}
    F_p2_F::Vector{Float64}
    F_p1_M::Vector{Float64}
    F_p2_M::Vector{Float64}
    F_p1_V::Vector{Float64}
    F_p2_V::Vector{Float64}
    F_p1_Ω::Vector{Float64}
    F_p2_Ω::Vector{Float64}
    F_F1_u::Vector{Float64}
    F_F2_u::Vector{Float64}
    F_F1_p::Vector{Float64}
    F_F2_p::Vector{Float64}
    F_F1_F::Vector{Float64}
    F_F2_F::Vector{Float64}
    F_F1_M::Vector{Float64}
    F_F2_M::Vector{Float64}
    F_M1_p::Vector{Float64}
    F_M2_p::Vector{Float64}
    F_M1_F::Vector{Float64}
    F_M2_F::Vector{Float64}
    F_M1_M::Vector{Float64}
    F_M2_M::Vector{Float64}
    F_V_u::Vector{Float64}
    F_V_p::Vector{Float64}
    F_V_V::Vector{Float64}
    F_Ω_p::Vector{Float64}
    F_Ω_Ω::Vector{Float64}
    F_u1_pdot::Vector{Float64}
    F_u2_pdot::Vector{Float64}
    F_u1_Vdot::Vector{Float64}
    F_u2_Vdot::Vector{Float64}
    F_u1_Ωdot::Vector{Float64}
    F_u2_Ωdot::Vector{Float64}
    F_p1_pdot::Vector{Float64}
    F_p2_pdot::Vector{Float64}
    F_p1_Vdot::Vector{Float64}
    F_p2_Vdot::Vector{Float64}
    F_p1_Ωdot::Vector{Float64}
    F_p2_Ωdot::Vector{Float64}
    F_V_udot::Vector{Float64}
    F_Ω_pdot::Vector{Float64}

    # Constructor
    function Jacobians()

        F_u1_p = zeros(3,3)
        F_u2_p = zeros(3,3)
        F_u1_F = zeros(3,3)
        F_u2_F = zeros(3,3)
        F_u1_V = zeros(3,3)
        F_u2_V = zeros(3,3)
        F_u1_Ω = zeros(3,3)
        F_u2_Ω = zeros(3,3)
        F_p1_p = zeros(3,3)
        F_p2_p = zeros(3,3)
        F_p1_F = zeros(3,3)
        F_p2_F = zeros(3,3)
        F_p1_M = zeros(3,3)
        F_p2_M = zeros(3,3)
        F_p1_V = zeros(3,3)
        F_p2_V = zeros(3,3)
        F_p1_Ω = zeros(3,3)
        F_p2_Ω = zeros(3,3)
        F_F1_u = zeros(3,3)
        F_F2_u = zeros(3,3)
        F_F1_p = zeros(3,3)
        F_F2_p = zeros(3,3)
        F_F1_F = zeros(3,3)
        F_F2_F = zeros(3,3)
        F_F1_M = zeros(3,3)
        F_F2_M = zeros(3,3)
        F_M1_p = zeros(3,3)
        F_M2_p = zeros(3,3)
        F_M1_F = zeros(3,3)
        F_M2_F = zeros(3,3)
        F_M1_M = zeros(3,3)
        F_M2_M = zeros(3,3)
        F_V_u = zeros(3,3)
        F_V_p = zeros(3,3)
        F_V_V = zeros(3,3)
        F_Ω_p = zeros(3,3)
        F_Ω_Ω = zeros(3,3)
        F_u1_pdot = zeros(3,3)
        F_u2_pdot = zeros(3,3)
        F_u1_Vdot = zeros(3,3)
        F_u2_Vdot = zeros(3,3)
        F_u1_Ωdot = zeros(3,3)
        F_u2_Ωdot = zeros(3,3)
        F_p1_pdot = zeros(3,3)
        F_p2_pdot = zeros(3,3)
        F_p1_Vdot = zeros(3,3)
        F_p2_Vdot = zeros(3,3)
        F_p1_Ωdot = zeros(3,3)
        F_p2_Ωdot = zeros(3,3)
        F_V_udot = zeros(3,3)
        F_Ω_pdot = zeros(3,3)

        return new(F_u1_p,F_u2_p,F_u1_F,F_u2_F,F_u1_V,F_u2_V,F_u1_Ω,F_u2_Ω,F_p1_p,F_p2_p,F_p1_F,F_p2_F,F_p1_M,F_p2_M,F_p1_V,F_p2_V,F_p1_Ω,F_p2_Ω,F_F1_u,F_F2_u,F_F1_p,F_F2_p,F_F1_F,F_F2_F,F_F1_M,F_F2_M,F_M1_p,F_M2_p,F_M1_F,F_M2_F,F_M1_M,F_M2_M,F_V_u,F_V_p,F_V_V,F_Ω_p,F_Ω_Ω,F_u1_pdot,F_u2_pdot,F_u1_Vdot,F_u2_Vdot,F_u1_Ωdot,F_u2_Ωdot,F_p1_pdot,F_p2_pdot,F_p1_Vdot,F_p2_Vdot,F_p1_Ωdot,F_p2_Ωdot,F_V_udot,F_Ω_pdot)
    end
end


#
# mutable struct Element <: BeamElement

#     Element composite type

#
mutable struct Element <: BeamElement

    # Parent beam
    parent::Beam
    # IDs
    localID::Int64
    globalID::Int64
    nodesLocalID::Vector{Int64}
    nodesGlobalID::Vector{Int64}
    # Point inertias attached to the element
    attachedPointInertias::Vector{PointInertia}
    # Geometric/material element variables
    Δℓ::Real
    x1::Real
    x1_norm::Float64
    k::Vector{<:Real}
    r::Vector{<:Real}
    R0::Matrix{Float64}
    R0T::Matrix{Float64}
    S::Matrix{Float64}
    I::Matrix{Float64}
    μ::Float64
    ηtilde::Matrix{Float64}
    S_11::Matrix{Float64}
    S_12::Matrix{Float64}
    S_21::Matrix{Float64}
    S_22::Matrix{Float64}
    I_11::Matrix{Float64}
    I_12::Matrix{Float64}
    I_21::Matrix{Float64}
    I_22::Matrix{Float64}
    # Geometric nodal variables
    x1_n1_norm::Float64
    x1_n2_norm::Float64
    x1_n1::Real
    x1_n2::Real
    r_n1::Vector{Float64}
    r_n2::Vector{Float64}
    R0_n1::Matrix{Float64}
    R0_n2::Matrix{Float64}
    R0T_n1::Matrix{Float64}
    R0T_n2::Matrix{Float64}
    # States, complementary states and their rates 
    states::ElementalStates{Float64}
    statesRates::ElementalStatesRates
    compStates::ComplementaryElementalStates{Float64}
    compStatesRates::ComplementaryElementalStatesRates
    # Nodal states 
    nodalStates::NodalStates{Float64}
    # Equivalent states' rates values at begin of time step
    udotEquiv::Vector{Float64}
    pdotEquiv::Vector{Float64}
    VdotEquiv::Vector{Float64}
    ΩdotEquiv::Vector{Float64}
    χdotEquiv::Vector{Float64}
    # Rotation variables
    R::Matrix{Float64}
    RR0::Matrix{Float64}
    RR0T::Matrix{Float64}
    RdotR0::Matrix{Float64}
    HT::Matrix{Float64}
    HTinv::Matrix{Float64}
    R_p1::Matrix{Float64}
    R_p2::Matrix{Float64}
    R_p3::Matrix{Float64}
    HT_p1::Matrix{Float64}
    HT_p2::Matrix{Float64}
    HT_p3::Matrix{Float64}
    HTinv_p1::Matrix{Float64}
    HTinv_p2::Matrix{Float64}
    HTinv_p3::Matrix{Float64}
    Rdot_p1::Matrix{Float64}
    Rdot_p2::Matrix{Float64}
    Rdot_p3::Matrix{Float64}
    # Velocities and accelerations of basis b
    v::Vector{Float64}
    ω::Vector{Float64}
    vdot::Vector{Float64}
    ωdot::Vector{Float64}
    # TF matrices for hinged nodes
    hingedNode1Mat::BitMatrix
    notHingedNode1Mat::BitMatrix
    notHingedNode2Mat::BitMatrix
    # Distributed loads' functions
    f_A_of_ζt::Function
    m_A_of_ζt::Function
    f_b_of_ζt::Function
    m_b_of_ζt::Function
    ff_A_of_ζt::Function
    mf_A_of_ζt::Function
    ff_b_of_ζt::Function
    mf_b_of_ζt::Function
    # Distributed loads' nodal values (over time)
    f_A::Array{Float64}
    m_A::Array{Float64}
    f_b::Array{Float64}
    m_b::Array{Float64}
    ff_A::Array{Float64}
    mf_A::Array{Float64}
    ff_b::Array{Float64}
    mf_b::Array{Float64}
    # Total distributed loads' nodal resultants
    f1::Vector{Float64}
    f2::Vector{Float64}
    m1::Vector{Float64}
    m2::Vector{Float64}
    # Gravitational loads' nodal resultants
    f_g::Vector{Float64}
    m_g::Vector{Float64}
    # Follower distributed loads' nodal resultants
    ff1_A::Vector{Float64}
    ff2_A::Vector{Float64}
    mf1_A::Vector{Float64}
    mf2_A::Vector{Float64}
    ff1_b::Vector{Float64}
    ff2_b::Vector{Float64}
    mf1_b::Vector{Float64}
    mf2_b::Vector{Float64}
    # Aerodynamic loads' nodal resultants
    f1_χ::Vector{Float64}
    f2_χ::Vector{Float64}
    m1_χ::Vector{Float64}
    m2_χ::Vector{Float64}
    # Derivatives w.r.t. extended rotation parameters
    f1_p::Matrix{Float64}
    f2_p::Matrix{Float64}
    m1_p::Matrix{Float64}
    m2_p::Matrix{Float64}
    # TFs for non-zero distributed loads
    hasDistributedDeadForcesBasisA::Bool
    hasDistributedDeadMomentsBasisA::Bool
    hasDistributedDeadForcesBasisb::Bool
    hasDistributedDeadMomentsBasisb::Bool
    hasDistributedFollowerForcesBasisA::Bool
    hasDistributedFollowerMomentsBasisA::Bool
    hasDistributedFollowerForcesBasisb::Bool
    hasDistributedFollowerMomentsBasisb::Bool
    # Indices on system of equations and related variables
    eqs_Fu1::Vector{Int64} 
    eqs_Fu2::Vector{Int64}
    eqs_Fp1::Vector{Int64}
    eqs_Fp2::Vector{Int64}
    eqs_FF1::Vector{Int64}
    eqs_FF2::Vector{Int64}
    eqs_FM1::Vector{Int64}
    eqs_FM2::Vector{Int64}
    eqs_FV::Vector{Int64}
    eqs_FΩ::Vector{Int64}
    eqs_Fχ::Vector{Int64}
    eqs_FF1_sep::Vector{Int64}
    eqs_FF2_sep::Vector{Int64}
    eqs_FM1_sep::Vector{Int64}
    eqs_FM2_sep::Vector{Int64}
    DOF_u::Vector{Int64}
    DOF_p::Vector{Int64}
    DOF_F::Vector{Int64}
    DOF_M::Vector{Int64}
    DOF_V::Vector{Int64}
    DOF_Ω::Vector{Int64}
    DOF_χ::Vector{Int64} 
    DOF_δ::Union{Vector{Int64},Int64}
    isSpecialNode1::Bool
    isSpecialNode2::Bool
    eqsNode1Set::Bool
    eqsNode2Set::Bool
    # Aerodynamic properties 
    aero::Union{Nothing,AeroProperties}
    # Relative rotation constraint
    rotationConstraint::Union{Nothing,RotationConstraint}
    # Finite element equation resultants and Jacobians
    resultants::Resultants
    jacobians::Jacobians

    # Constructor
    function Element(parent::Beam) 

        # ID on beam (local)
        localID = length(parent.elements) + 1

        # ID on assembly (global - updated later on model definition)
        globalID = localID

        # Nodes' IDs on beam (local)
        nodesLocalID = [localID, localID+1]

        # Nodes' IDs on assembly (global - updated later on model definition)
        nodesGlobalID = nodesLocalID

        # Initialize point inertias attached to the element
        attachedPointInertias = Vector{PointInertia}()

        # Nodal arclength coordinates normalized by parent beam length
        x1_n1_norm,x1_n2_norm = parent.normalizedNodalPositions[nodesLocalID]

        # Nodal arclength positions
        x1_n1,x1_n2 = round_off!(parent.length * parent.normalizedNodalPositions[nodesLocalID])

        # Length
        Δℓ = x1_n2-x1_n1

        # Midpoint arclength position
        x1 = (x1_n1+x1_n2)/2

        # Midpoint arclength position normalized by parent beam length
        x1_norm = x1/parent.length

        # Initial curvature vector
        k = parent.k

        # Get nodal coordinates and add to parent beam
        r_n1 = position_vector_from_curvature(parent.R0,k,x1_n1)
        r_n2 = position_vector_from_curvature(parent.R0,k,x1_n2)

        # Nodal rotation tensors
        R0_n1 = rotation_tensor_from_curvature(parent.R0,k,x1_n1)
        R0_n2 = rotation_tensor_from_curvature(parent.R0,k,x1_n2)
        R0T_n1 = Matrix(R0_n1')
        R0T_n2 = Matrix(R0_n2')

        # Midpoint coordinates
        r = position_vector_from_curvature(parent.R0,k,x1)

        # Midpoint rotation tensor from basis A to basis b and its transpose
        R0 = rotation_tensor_from_curvature(parent.R0,k,x1)
        R0T = Matrix(R0')

        # Sectional compliance (force-strain) matrix and submatrices
        if length(parent.C) == 1
            S = parent.C[1]^-1
        else
            S = parent.C[localID]^-1
        end
        S_11 = S[1:3,1:3]
        S_12 = S[1:3,4:6]
        S_21 = S[4:6,1:3]
        S_22 = S[4:6,4:6]
        
        # Sectional inertia (momentum-velocity) matrix (includes contributions from point inertias)
        if length(parent.I) == 1
            I = parent.I[1]
        else
            I = parent.I[localID]
        end

        # Loop point inertias attached to the beam
        if !isnothing(parent.pointInertias)
            for pointInertia in parent.pointInertias
                # Check if is linked to current element
                if localID == pointInertia.elementID 
                    @unpack mass,η,inertiaMatrix = pointInertia
                    # Increment element's sectional inertia matrix
                    I += 1/Δℓ * [mass*I3 -mass*tilde(η);
                        mass*tilde(η) inertiaMatrix-mass*tilde(η)^2]
                    # Add to element's list of attached point inertias
                    push!(attachedPointInertias,pointInertia)    
                end
            end 
        end

        # Sectional inertia submatrices (includes contributions from point inertias)
        I_11 = I[1:3,1:3]
        I_12 = I[1:3,4:6]
        I_21 = I[4:6,1:3]
        I_22 = I[4:6,4:6]

        # Mass per unit length and sectional centroid position's skew-symmetric matrix with respect to reference point, resolved in basis b
        μ = I[1,1]
        ηtilde = μ > 0 ? I_21/μ : zeros(3,3)
        
        # Initialize states, complementary states and states' rates 
        states = ElementalStates{Float64}()
        statesRates = ElementalStatesRates()
        compStates = ComplementaryElementalStates{Float64}()
        compStatesRates = ComplementaryElementalStatesRates()
        nodalStates = NodalStates{Float64}()

        # Initialize equivalent states' rates at begin of time step
        udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv = zeros(3),zeros(3),zeros(3),zeros(3)
        χdotEquiv = zeros(0)

        # Initialize rotation tensors
        R,RR0,RR0T,HT,HTinv,RdotR0,R_p1,R_p2,R_p3,HT_p1,HT_p2,HT_p3,HTinv_p1,HTinv_p2,HTinv_p3,Rdot_p1,Rdot_p2,Rdot_p3 = I3,R0,R0T,R0,I3,I3,zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3)

        # Initialize element velocities and accelerations in basis b
        v,ω,vdot,ωdot = zeros(3),zeros(3),zeros(3),zeros(3)

        # Get matrices of hinged and not hinged nodes' DoF
        hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat = get_hinged_nodes_matrices(parent,nodesLocalID)

        # Get distributed loads as a function of the local coordinate
        f_A_of_ζt,m_A_of_ζt,f_b_of_ζt,m_b_of_ζt,ff_A_of_ζt,mf_A_of_ζt,ff_b_of_ζt,mf_b_of_ζt = get_element_distributed_loads(parent,Δℓ,x1_n1)

        # Initialize nodal values of distributed loads 
        f_A,m_A,f_b,m_b,ff_A,mf_A,ff_b,mf_b = zeros(2,3,1),zeros(2,3,1),zeros(2,3,1),zeros(2,3,1),zeros(2,3,1),zeros(2,3,1),zeros(2,3,1),zeros(2,3,1)

        # Initialize nodal resultants of distributed loads, resolved in basis A
        f1,f2,m1,m2,f_g,m_g,ff1_A,ff2_A,mf1_A,mf2_A,ff1_b,ff2_b,mf1_b,mf2_b,f1_χ,f2_χ,m1_χ,m2_χ = zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3)
    
        # Initialize derivatives of total nodal resultants of distributed loads w.r.t. extended rotation parameters
        f1_p,f2_p,m1_p,m2_p = zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3)

        # Initialize TFs for non-zero nodal values of distributed loads 
        hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb = false,false,false,false,false,false,false,false

        # Initialize system indices
        eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,eqs_Fχ,eqs_FF1_sep,eqs_FF2_sep,eqs_FM1_sep,eqs_FM2_sep,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,DOF_χ,DOF_δ = Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}()

        isSpecialNode1,isSpecialNode2,eqsNode1Set,eqsNode2Set = false,false,false,false

        # Initialize aerodynamic properties 
        aero = !isnothing(parent.aeroSurface) ? AeroProperties(parent.aeroSurface,R0,x1,x1_norm,x1_n1_norm,x1_n2_norm) : nothing

        # Update aerodynamic states and rates vectors to correct size
        if !isnothing(aero)
            states.χ = zeros(aero.nTotalAeroStates)
            if typeof(aero.solver) == BLi
                states.χ[aero.pitchPlungeStatesRange[4:6]] .= 1.0
            elseif typeof(aero.solver) == BLo
                states.χ[10:11] .= 1.0
            end
            statesRates.χdot = zeros(aero.nTotalAeroStates)
            χdotEquiv = zeros(aero.nTotalAeroStates)
        end

        # Initialize rotation constraint on slave element
        rotationConstraint = nothing

        # Initialize equation resultants and Jacobians
        resultants = Resultants(aero.nTotalAeroStates)
        jacobians = Jacobians()
        
        # Create element
        self = new(parent,localID,globalID,nodesLocalID,nodesGlobalID,attachedPointInertias,Δℓ,x1,x1_norm,k,r,R0,R0T,S,I,μ,ηtilde,S_11,S_12,S_21,S_22,I_11,I_12,I_21,I_22,x1_n1_norm,x1_n2_norm,x1_n1,x1_n2,r_n1,r_n2,R0_n1,R0_n2,R0T_n1,R0T_n2,states,statesRates,compStates,compStatesRates,nodalStates,udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv,χdotEquiv,R,RR0,RR0T,RdotR0,HT,HTinv,R_p1,R_p2,R_p3,HT_p1,HT_p2,HT_p3,HTinv_p1,HTinv_p2,HTinv_p3,Rdot_p1,Rdot_p2,Rdot_p3,v,ω,vdot,ωdot,hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat,f_A_of_ζt,m_A_of_ζt,f_b_of_ζt,m_b_of_ζt,ff_A_of_ζt,mf_A_of_ζt,ff_b_of_ζt,mf_b_of_ζt,f_A,m_A,f_b,m_b,ff_A,mf_A,ff_b,mf_b,f1,f2,m1,m2,f_g,m_g,ff1_A,ff2_A,mf1_A,mf2_A,ff1_b,ff2_b,mf1_b,mf2_b,f1_χ,f2_χ,m1_χ,m2_χ,f1_p,f2_p,m1_p,m2_p,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb,eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,eqs_Fχ,eqs_FF1_sep,eqs_FF2_sep,eqs_FM1_sep,eqs_FM2_sep,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,DOF_χ,DOF_δ,isSpecialNode1,isSpecialNode2,eqsNode1Set,eqsNode2Set,aero,rotationConstraint,resultants,jacobians)

        # Add element to parent beam
        push!(parent.elements, self)

        return self
    end
end


# Gets the TF matrices resulting from hinged and not hinged nodal DoFs times the identity matrix
function get_hinged_nodes_matrices(parent::Beam,nodesLocalID::Vector{Int64})

    @unpack hingedNodes,hingedNodesDoF = parent

    I3bit = BitMatrix(Matrix(true*LinearAlgebra.I,3,3))

    # Initialize outputs
    hingedNode1Mat,hingedNode2Mat,notHingedNode1Mat = falses(3,3),falses(3,3),I3bit

    # Initialize TF of hinged DoFs on element's nodes
    hingedNode1DoF,hingedNode2DoF = falses(3),falses(3)

    # Loop over beam's hinged nodes
    for (hingedNode,hingedNodeDoF) in zip(hingedNodes,hingedNodesDoF)
        # Skip if the element does not contain the hinged node
        if !(hingedNode in nodesLocalID)
            continue
        end
        # TF of hinged DoFs on element's first node
        hingedNode1DoF = nodesLocalID[1] == hingedNode ? hingedNodeDoF : falses(3)
        # TF of hinged DoFs on element's second node
        hingedNode2DoF = nodesLocalID[2] == hingedNode ? hingedNodeDoF : falses(3)
        break
    end

    # Hinged node TF matrices
    hingedNode1Mat = (hingedNode1DoF*hingedNode1DoF') .* I3bit
    hingedNode2Mat = (hingedNode2DoF*hingedNode2DoF') .* I3bit
    notHingedNode1Mat = .!hingedNode1Mat .* I3bit
    notHingedNode2Mat = .!hingedNode2Mat .* I3bit

    return hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat
end


# Gets the distributed loads in the element's local coordinate
function get_element_distributed_loads(parent::Beam,Δℓ::Real,x1_n1::Real)

    @unpack f_A_of_x1t,m_A_of_x1t,f_b_of_x1t,m_b_of_x1t,ff_A_of_x1t,mf_A_of_x1t,ff_b_of_x1t,mf_b_of_x1t = parent

    # Function for arclength in terms of local coordinate ζ
    x1 = ζ -> x1_n1 + Δℓ*ζ

    # Dead loads in basis A
    f_A_of_ζt = (ζ,t) -> f_A_of_x1t(x1(ζ),t)
    m_A_of_ζt = (ζ,t) -> m_A_of_x1t(x1(ζ),t)

    # Dead loads in basis b
    f_b_of_ζt = (ζ,t) -> f_b_of_x1t(x1(ζ),t)
    m_b_of_ζt = (ζ,t) -> m_b_of_x1t(x1(ζ),t)

    # Follower loads in basis A
    ff_A_of_ζt = (ζ,t) -> ff_A_of_x1t(x1(ζ),t)
    mf_A_of_ζt = (ζ,t) -> mf_A_of_x1t(x1(ζ),t)

    # Follower loads in basis b
    ff_b_of_ζt = (ζ,t) -> ff_b_of_x1t(x1(ζ),t)
    mf_b_of_ζt = (ζ,t) -> mf_b_of_x1t(x1(ζ),t)

    return f_A_of_ζt,m_A_of_ζt,f_b_of_ζt,m_b_of_ζt,ff_A_of_ζt,mf_A_of_ζt,ff_b_of_ζt,mf_b_of_ζt
end


# Update the distributed loads in the element's local coordinate
function update_element_distributed_loads!(element::Element,loadType::String,loadFun::Function)

    @unpack x1_n1,Δℓ = element

    # Function for arclength in terms of local coordinate ζ
    x1(ζ) = x1_n1 + Δℓ*ζ

    # Update loads
    if loadType == "f_A_of_x1t"
        element.f_A_of_ζt = (ζ,t) -> loadFun(x1(ζ),t)
    elseif loadType == "m_A_of_x1t"
        element.m_A_of_ζt = (ζ,t) -> loadFun(x1(ζ),t)
    elseif loadType == "f_b_of_x1t"
        element.f_b_of_ζt = (ζ,t) -> loadFun(x1(ζ),t)
    elseif loadType == "m_b_of_x1t"
        element.m_b_of_ζt = (ζ,t) -> loadFun(x1(ζ),t)
    elseif loadType == "ff_A_of_x1t"
        element.ff_A_of_ζt = (ζ,t) -> loadFun(x1(ζ),t)
    elseif loadType == "mf_A_of_x1t"
        element.mf_A_of_ζt = (ζ,t) -> loadFun(x1(ζ),t)
    elseif loadType == "ff_b_of_x1t"
        element.ff_b_of_ζt = (ζ,t) -> loadFun(x1(ζ),t)
    elseif loadType == "mf_b_of_x1t"
        element.mf_b_of_ζt = (ζ,t) -> loadFun(x1(ζ),t)
    end

end


# Adds the point inertia's matrix to the element's sectional inertia matrix
function add_point_inertia_to_element!(element::Element,pointInertia::PointInertia)

    @unpack I,Δℓ,attachedPointInertias = element
    @unpack η,mass,inertiaMatrix = pointInertia

    # Add to inertia matrix 
    I += 1/Δℓ * [mass*I3 -mass*tilde(η);
    mass*tilde(η) inertiaMatrix-mass*tilde(η)^2]

    # Update submatrices
    I_11 = I[1:3,1:3]
    I_12 = I[1:3,4:6]
    I_21 = I[4:6,1:3]
    I_22 = I[4:6,4:6]

    # Update mass per unit length and skew-symmetric matrix of sectional position
    μ = I[1,1]
    ηtilde = μ > 0 ? I_21/μ : zeros(3,3)

    # Add to element's list of attached point inertias
    push!(attachedPointInertias,pointInertia)

    @pack! element = I,I_11,I_12,I_21,I_22,μ,ηtilde,attachedPointInertias

end
