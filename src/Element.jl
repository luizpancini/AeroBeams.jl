"""
mutable struct ElementalStates

    ElementalStates composite type

# Fields

"""
mutable struct ElementalStates{T<:Union{Float64,ComplexF64}}
    
    # Fields
    u::Vector{T}
    p::Vector{T}
    F::Vector{T}
    M::Vector{T}
    V::Vector{T}
    Ω::Vector{T}

    # Constructor
    function ElementalStates{T}() where T<:Union{Float64,ComplexF64}

        u,p,F,M,V,Ω = zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3)

        return new{T}(u,p,F,M,V,Ω)
    end
end


"""
mutable struct ComplementaryElementalStates

    ComplementaryElementalStates composite type

# Fields

"""
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
    
"""
mutable struct ElementalStatesRates

    ElementalStatesRates composite type

# Fields

"""
mutable struct ElementalStatesRates
    
    # Fields
    udot::Vector{Float64}
    pdot::Vector{Float64}
    Vdot::Vector{Float64}
    Ωdot::Vector{Float64}
    uddot::Vector{Float64}
    pddot::Vector{Float64}

    # Constructor
    function ElementalStatesRates() 

        udot,pdot,Vdot,Ωdot,uddot,pddot = zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3)

        return new(udot,pdot,Vdot,Ωdot,uddot,pddot)
    end
end


"""
mutable struct ComplementaryElementalStatesRates

    ComplementaryElementalStatesRates composite type

# Fields

"""
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


"""
mutable struct NodalStates

    NodalStates composite type

# Fields

"""
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


"""
mutable struct Element <: BeamElement

    Element composite type

# Fields
- parent::Beam
- localID::Int64
# Notes
 - Finite elements belong to a Beam
"""
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
    Δℓ::Float64
    x1::Float64
    k::Vector{Float64}
    r::Vector{Float64}
    R0::Matrix{Float64}
    R0T::Matrix{Float64}
    R_cs::Matrix{Float64}
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
    x1_n1::Float64
    x1_n2::Float64
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
    uddotEquiv::Vector{Float64}
    pddotEquiv::Vector{Float64}
    # Rotation tensor functions
    R::Matrix{Float64}
    RR0::Matrix{Float64}
    RR0T::Matrix{Float64}
    RRwR0::Matrix{Float64}
    RdotR0::Matrix{Float64}
    HT::Matrix{Float64}
    HTinv::Matrix{Float64}
    # Velocities and accelerations in basis b
    v_b::Vector{Float64}
    ω_b::Vector{Float64}
    vdot_b::Vector{Float64}
    ωdot_b::Vector{Float64}
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
    # Current values of nodal resultants
    f1::Vector{Float64}
    f2::Vector{Float64}
    m1::Vector{Float64}
    m2::Vector{Float64}
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
    isSpecialNode1::Bool
    isSpecialNode2::Bool
    eqsNode1Set::Bool
    eqsNode2Set::Bool

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

        # Nodal arclength positions
        x1_n1,x1_n2 = parent.length .* parent.normalizedNodalPositions[nodesLocalID]

        # Length
        Δℓ = x1_n2-x1_n1

        # Midpoint arclength position
        x1 = (x1_n1+x1_n2)/2

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

        # Midpoint rotation tensor and its transpose
        R0 = rotation_tensor_from_curvature(parent.R0,k,x1)
        R0T = Matrix(R0')

        # Midpoint cross-section rotation tensor
        R_cs = rotation_tensor_from_curvature(parent.R_cs,k,x1)

        # Compound cross-section rotation tensor
        R_cs_comp = [R_cs' zeros(3,3);
                     zeros(3,3) R_cs']

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
        for pointInertia in parent.pointInertias
            # Check if is linked to current element
            if localID == pointInertia.elementLocalID 
                @unpack mass,ηtilde,inertiaMatrix = pointInertia
                # Increment element's sectional inertia matrix
                I += 1/Δℓ * [mass*I3 -mass*ηtilde;
                    mass*ηtilde inertiaMatrix-mass*ηtilde^2]
                # Add to element's list of attached point inertias
                push!(attachedPointInertias,pointInertia)    
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
        udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv,uddotEquiv,pddotEquiv = zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3)

        # Initialize rotation tensors
        R,RR0,RR0T,RRwR0,HT,HTinv,RdotR0 = I3,R0,R0T,R0,I3,I3,zeros(3,3)

        # Initialize element velocities and accelerations in basis b
        v_b,ω_b,vdot_b,ωdot_b = zeros(3),zeros(3),zeros(3),zeros(3)

        # Get matrices of hinged and not hinged nodes' DoF
        hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat = get_hinged_nodes_matrices(parent,nodesLocalID)

        # Get distributed loads as a function of the local coordinate
        f_A_of_ζt,m_A_of_ζt,f_b_of_ζt,m_b_of_ζt,ff_A_of_ζt,mf_A_of_ζt,ff_b_of_ζt,mf_b_of_ζt = get_element_distributed_loads(parent,Δℓ,x1_n1)

        # Initialize nodal values of distributed loads 
        f_A,m_A,f_b,m_b,ff_A,mf_A,ff_b,mf_b = zeros(2,3,1),zeros(2,3,1),zeros(2,3,1),zeros(2,3,1),zeros(2,3,1),zeros(2,3,1),zeros(2,3,1),zeros(2,3,1)

        # Initialize resultants of nodal generalized forces, resolved in basis A
        f1,f2,m1,m2 = zeros(3),zeros(3),zeros(3),zeros(3)

        # Initialize TFs for non-zero nodal values of distributed loads 
        hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb = false,false,false,false,false,false,false,false

        # Initialize system indices
        eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,eqs_FF1_sep,eqs_FF2_sep,eqs_FM1_sep,eqs_FM2_sep,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω = Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}(),Vector{Int64}()

        isSpecialNode1,isSpecialNode2,eqsNode1Set,eqsNode2Set = false,false,false,false
        
        # Create element
        self = new(parent,localID,globalID,nodesLocalID,nodesGlobalID,attachedPointInertias,Δℓ,x1,k,r,R0,R0T,R_cs,S,I,μ,ηtilde,S_11,S_12,S_21,S_22,I_11,I_12,I_21,I_22,x1_n1,x1_n2,r_n1,r_n2,R0_n1,R0_n2,R0T_n1,R0T_n2,states,statesRates,compStates,compStatesRates,nodalStates,udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv,uddotEquiv,pddotEquiv,R,RR0,RR0T,RRwR0,RdotR0,HT,HTinv,v_b,ω_b,vdot_b,ωdot_b,hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat,f_A_of_ζt,m_A_of_ζt,f_b_of_ζt,m_b_of_ζt,ff_A_of_ζt,mf_A_of_ζt,ff_b_of_ζt,mf_b_of_ζt,f_A,m_A,f_b,m_b,ff_A,mf_A,ff_b,mf_b,f1,f2,m1,m2,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb,eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,eqs_FF1_sep,eqs_FF2_sep,eqs_FM1_sep,eqs_FM2_sep,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,isSpecialNode1,isSpecialNode2,eqsNode1Set,eqsNode2Set)

        # Add element to parent beam
        push!(parent.elements, self)

        return self
    end
end


"""
get_hinged_nodes_matrices(parent::Beam,nodesLocalID::Vector{Int64})

Gets the TF matrices resulting from hinged and not hinged nodal DoFs times the identity matrix

# Arguments
- parent::Beam
- nodesLocalID::Vector{Int64}
"""
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


"""
get_element_distributed_loads(parent::Beam,Δℓ::Float64,x1_n1::Float64)

Gets the distributed loads in the element's local coordinate

# Arguments
- parent::Beam
- Δℓ::Float64
- x1_n1::Float64
"""
function get_element_distributed_loads(parent::Beam,Δℓ::Float64,x1_n1::Float64)

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


"""
update_element_distributed_loads!(element::Element,loadType::String,loadFun::Function)

Update the distributed loads in the element's local coordinate

# Arguments
- element::Element
- loadType::String
- loadFun::Function
"""
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


"""
add_point_inertia_to_element!(element::Element,pointInertia::PointInertia)

Adds the point inertia's matrix to the element's sectional inertia matrix

# Arguments
- element::Element
- pointInertia::PointInertia
"""
function add_point_inertia_to_element!(element::Element,pointInertia::PointInertia)

    @unpack I,Δℓ = element
    @unpack mass,ηtilde,inertiaMatrix = pointInertia

    I += 1/Δℓ * [mass*I3 -mass*ηtilde;
    mass*ηtilde inertiaMatrix-mass*ηtilde^2]

    @pack! element = I

end
