## MODEL COMPOSITE TYPES
# ------------------------------------------------------------------------------
# KinematicsBasisAState composite type
# Fields
# - `u_A`: current displacement of basis A relative to basis I
# - `v_A`: current velocity of basis A relative to basis I
# - `ω_A`: current angular velocity of basis A relative to basis I
# - `vdot_A`: current acceleration of basis A relative to basis I
# - `ωdot_A`: current angular acceleration of basis A relative to basis I
# - `R_A`: current rotation tensor from basis I to basis A
# - `R_AT`: current rotation tensor from basis A to basis I
#
struct KinematicsBasisAState

    u_A::MVector{3,Float64}
    v_A::MVector{3,Float64}
    ω_A::MVector{3,Float64}
    vdot_A::MVector{3,Float64}
    ωdot_A::MVector{3,Float64}
    R_A::MMatrix{3,3,Float64,9}
    R_AT::MMatrix{3,3,Float64,9}

end
KinematicsBasisAState(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A) = KinematicsBasisAState(
    MVector{3,Float64}(u_A(0)),
    MVector{3,Float64}(v_A(0)),
    MVector{3,Float64}(ω_A(0)),
    MVector{3,Float64}(vdot_A(0)),
    MVector{3,Float64}(ωdot_A(0)),
    MMatrix{3,3,Float64,9}(rotation_tensor_WM(p_A0)),
    MMatrix{3,3,Float64,9}(rotation_tensor_WM(p_A0)'),
)


# KinematicsBasisA composite type
# Fields
# - `p_A0`: initial Wiener-Milenković rotation parameters that bring basis I to basis A
# - `u_A`: displacement of basis A relative to basis I
# - `v_A`: velocity of basis A relative to basis I
# - `ω_A`: angular velocity of basis A relative to basis I
# - `vdot_A`: acceleration of basis A relative to basis I
# - `ωdot_A`: angular acceleration of basis A relative to basis I
# - `u_AVec`: time vector of displacement of basis A relative to basis I
# - `v_AVec`: time vector of velocity of basis A relative to basis I
# - `ω_AVec`: time vector of angular velocity of basis A relative to basis I
# - `vdot_AVec`: time vector of acceleration of basis A relative to basis I
# - `ωdot_AVec`: time vector of angular acceleration of basis A relative to basis I
# - `R_AVec`: time vector of rotation tensor from basis I to basis A
# - `state`: state of the kinematics of basis A
#
mutable struct KinematicsBasisA

    p_A0::SVector{3,Float64}
    u_A::Function
    v_A::Function
    ω_A::Function
    vdot_A::Function
    ωdot_A::Function

    u_AVec::Vector{SVector{3,Float64}}
    v_AVec::Vector{SVector{3,Float64}}
    ω_AVec::Vector{SVector{3,Float64}}
    vdot_AVec::Vector{SVector{3,Float64}}
    ωdot_AVec::Vector{SVector{3,Float64}}
    R_AVec::Vector{SMatrix{3,3,Float64,9}}

    state::KinematicsBasisAState

end
KinematicsBasisA(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A) = 
    KinematicsBasisA(
        SVector{3,Float64}(p_A0), u_A, v_A, ω_A, vdot_A, ωdot_A,
        [SVector{3,Float64}(u_A(0))], [SVector{3,Float64}(v_A(0))], [SVector{3,Float64}(ω_A(0))],
        [SVector{3,Float64}(vdot_A(0))], [SVector{3,Float64}(ωdot_A(0))], [SMatrix{3,3,Float64,9}(rotation_tensor_WM(p_A0))],
        KinematicsBasisAState(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A),
    )


# InitialInertiaProperties composite type
# Fields
# - `mass`: total mass of the model
# - `centerOfMass`: center of mass of the model, resolved in basis A
# - `inertiaMatrix`: inertia matrix of the model about its center of mass, resolved in basis A
#
struct InitialInertiaProperties

    mass::Float64
    centerOfMass::SVector{3,Float64}
    inertiaMatrix::SMatrix{3,3,Float64,9}

end


# Model assembly data composite type
# Fields
# - `elements`: elements in the model (beams' elements)
# - `elementGroups`: the same elements, in the same order, split into a tuple of concretely-typed vectors (consecutive elements of the same type); allows type-stable loops over the elements
# - `nElementsTotal`: total number of elements in the model
# - `nNodesTotal`: total number of nodes in the model
# - `elementNodes`: global IDs of the nodes of each element in the model
# - `r_n`: nodal coordinates of the nodes of each element in the model, relative to basis A
# - `specialNodes`: special nodes in the model
# - `specialNodesGlobalIDs`: global IDs of the special nodes in the model
# - `BCedNodes`: global IDs of the BCed nodes in the model
# - `springedNodes`: global IDs of the springed nodes in the model
# - `springsOnNodes`: springs attached to the springed nodes in the model
# - `systemOrder`: order of the system (number of DOFs)
# - `nTrimVariables`: total number of trim variables in the model
# - `nTrimLoads`: number of trim loads in the model
# - `nTrimδ`: number of trim flap deflections in the model
# - `rangeNominal`: range of indices corresponding to the regular system in the system state vector
# - `rangeAeroStates`: range of indices corresponding all the aerodynamic states in the model
# - `rangeTrimLoads`: range of indices corresponding to the trim loads in the system state vector
# - `rangeTrimδ`: range of indices corresponding to the trim flap deflections in the system state vector
# - `forceScaling`: scaling factor for the forces in the model
# - `invForceScaling`: inverse of the above
# - `initialInertiaProperties`: initial inertia properties of the model
#
struct ModelAssemblyData

    elements::Vector{Element}
    elementGroups::Tuple
    nElementsTotal::Int
    nNodesTotal::Int
    elementNodes::Vector{SVector{2,Int}}
    r_n::Vector{SVector{3,Float64}}

    forceScaling::Float64
    invForceScaling::Float64

    BCedNodes::Vector{Int}
    specialNodes::Vector{SpecialNode}
    specialNodesGlobalIDs::Vector{Int}

    systemOrder::Int
    nTrimVariables::Int
    nTrimLoads::Int
    nTrimδ::Int
    rangeNominal::UnitRange{Int}
    rangeAeroStates::Vector{Int}
    rangeTrimLoads::UnitRange{Int}
    rangeTrimδ::UnitRange{Int}

    initialInertiaProperties::InitialInertiaProperties

end


# InitialModelState composite type
# Fields
# - `elementsStates`: states of the elements in the model
# - `specialNodesStates`: states of the special nodes in the model
#
struct InitialModelState

    elementsStates::Vector{ElementState}
    specialNodesStates::Vector{SpecialNodeState}

end


# Model composite type
# Fields
# - `name`: name of the model
# - `units`: units system
# - `firstNodePosition`: position vector of the first node of the first beam, resolved in the inertial frame I
# - `gravityVector`: gravity vector, resolved in the inertial frame I
# - `beams`: beams in the assembly
# - `connections`: connections between beams and nodes in the assembly
# - `BCs`: boundary condtions
# - `kinematicsBasisA`: kinematics basis A
# - `trimLoadsLinks`: trim loads links
# - `flapLinks`: flap links
# - `hingeAxisConstraints`: hinge axis constraints
# - `atmosphere`: atmosphere
# - `gust`: gust
# - `assemblyData`: assembly data for the model
#
@with_kw mutable struct Model

    # Definition fields (inputs for Model creation)
    name::String
    units::UnitsSystem
    firstNodePosition::SVector{3,Float64}
    gravityVector::SVector{3,Float64}
    beams::Vector{Beam}
    connections::Vector{Tuple{Pair{Beam,Beam},Pair{Int,Int}}}
    BCs::Vector{BC}
    kinematicsBasisA::KinematicsBasisA
    trimLoadsLinks::Vector{TrimLoadsLink}
    flapLinks::Vector{FlapLink}
    hingeAxisConstraints::Vector{<:HingeAxisConstraint}
    atmosphere::Atmosphere
    gust::Union{Nothing,Gust}

    # Assembly field (deduced from inputs)
    assemblyData::ModelAssemblyData

    # Initial state field
    initialState::InitialModelState

end
export Model


## PUBLIC API
# ------------------------------------------------------------------------------
"""
    create_Model(; kwargs...)

Creates a model

# Keyword arguments
- `name::String`: name of the model
- `units::UnitsSystem`: units system
- `firstNodePosition::AbstractVector{<:Real}`: position vector from basis A to first node of the first beam, resolved in basis A
- `gravityVector::AbstractVector{<:Real}`: gravity vector, resolved in the inertial frame I
- `beams::Vector{Beam}`: beams in the assembly
- `connections::Vector{<:Tuple{<:Pair{Beam,Beam},<:Pair{Int,Int}}}`: beam/node connection pairs (`(beamThis => beamOther, nodeThis => nodeOther)`)
- `BCs::Vector{BC}`: boundary condtions
- `p_A0::AbstractVector{<:Real}`: initial Wiener-Milenković rotation parameters that bring basis I to basis A
- `u_A::Union{Nothing,AbstractVector{<:Real},<:Function}`: displacement of basis A relative to basis I
- `v_A::Union{Nothing,AbstractVector{<:Real},<:Function}`: velocity of basis A relative to basis I
- `ω_A::Union{Nothing,AbstractVector{<:Real},<:Function}`: angular velocity of basis A relative to basis I
- `vdot_A::Union{Nothing,AbstractVector{<:Real},<:Function}`: acceleration of basis A relative to basis I
- `ωdot_A::Union{Nothing,AbstractVector{<:Real},<:Function}`: angular acceleration of basis A relative to basis I
- `trimLoadsLinks::Vector{TrimLoadsLink}`: links between trim loads
- `flapLinks::Vector{FlapLink}`: links between flapped surfaces
- `hingeAxisConstraints::Vector{<:HingeAxisConstraint}`: hinge axis constraints
- `atmosphere::Atmosphere`: atmosphere
- `gust::Union{Nothing,Gust}`: gust
"""
function create_Model(; 
    name::String="",
    units::UnitsSystem=create_UnitsSystem(),
    firstNodePosition::AbstractVector{<:Real}=zeros(3),
    gravityVector::AbstractVector{<:Real}=zeros(3),
    beams::Vector{Beam},
    connections::Vector{<:Tuple{<:Pair{Beam,Beam},<:Pair{Int,Int}}}=Tuple{Pair{Beam,Beam},Pair{Int,Int}}[],
    BCs::Vector{BC}=Vector{BC}(),
    p_A0::AbstractVector{<:Real}=zeros(3),
    u_A::Union{Nothing,AbstractVector{<:Real},<:Function}=nothing,
    v_A::Union{Nothing,AbstractVector{<:Real},<:Function}=nothing,
    ω_A::Union{Nothing,AbstractVector{<:Real},<:Function}=nothing,
    vdot_A::Union{Nothing,AbstractVector{<:Real},<:Function}=nothing,
    ωdot_A::Union{Nothing,AbstractVector{<:Real},<:Function}=nothing,
    trimLoadsLinks::Vector{TrimLoadsLink}=Vector{TrimLoadsLink}(),
    flapLinks::Vector{FlapLink}=Vector{FlapLink}(),
    hingeAxisConstraints::Vector{<:HingeAxisConstraint}=HingeAxisConstraint[],
    atmosphere::Atmosphere=standard_atmosphere(0),
    gust::Union{Nothing,Gust}=nothing
    )

    # Validate inputs
    firstNodePosition,gravityVector,connections,p_A0,uAIsInput,vAIsInput,ωAIsInput,vdotAIsInput,ωdotAIsInput = validate_model(beams,firstNodePosition,gravityVector,connections,p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A)

    # Compute motion of basis A and set kinematics object
    kinematicsBasisA = get_kinematics_basis_A(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A,uAIsInput,vAIsInput,ωAIsInput,vdotAIsInput,ωdotAIsInput)

    # Get assembly data and update BCs
    assemblyData,BCs = get_model_assembly_data!(beams,connections,firstNodePosition,trimLoadsLinks,flapLinks,hingeAxisConstraints,BCs,gust)

    # Update initial states of elements and special nodes
    initialState = update_initial_states!(beams,assemblyData.elements,assemblyData.specialNodes,kinematicsBasisA,atmosphere)
    
    return Model(
        name=name,
        units=units,
        firstNodePosition=firstNodePosition,
        gravityVector=gravityVector,
        beams=beams,
        connections=connections,
        BCs=BCs,
        kinematicsBasisA=kinematicsBasisA,
        trimLoadsLinks=trimLoadsLinks,
        flapLinks=flapLinks,
        hingeAxisConstraints=hingeAxisConstraints,
        atmosphere=atmosphere,
        gust=gust,
        assemblyData=assemblyData,
        initialState=initialState,
        )

end
export create_Model


"""
    update_initial_rotation_of_basis_A!(model; p_A0)

Updates the initial rotation of basis A (and related kinematics) in the model.

# Arguments
- `model::Model`: reference model

# Keyword Arguments
- `p_A0::AbstractVector{<:Real}`: initial rotation of basis A
"""
function update_initial_rotation_of_basis_A!(model::Model; p_A0::AbstractVector{<:Real}) 

    # Validate input and set correct type
    @assert length(p_A0) == 3 "p_A0 must be a 3-element vector"
    p_A0 = SVector{3,Float64}(p_A0)

    # Unpack kinematics of basis A
    @unpack u_A,v_A,ω_A,vdot_A,ωdot_A = model.kinematicsBasisA

    # Get new kinematics of basis A
    kinematicsBasisA = KinematicsBasisA(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A)

    # Update initial states
    initialState = update_initial_states!(model.beams,model.assemblyData.elements,model.assemblyData.specialNodes,kinematicsBasisA,model.atmosphere)

    # Set on model
    setfield!(model, :kinematicsBasisA, kinematicsBasisA)
    setfield!(model, :initialState, initialState)

end
export update_initial_rotation_of_basis_A!


"""
    update_u_A!(model; u_A)

Updates the displacement of basis A (and related kinematics) in the model.

# Arguments
- `model::Model`: reference model

# Keyword Arguments
- `u_A::Union{AbstractVector{<:Real},<:Function}`: displacement of basis A
"""
function update_u_A!(model::Model; u_A::Union{AbstractVector{<:Real},<:Function})

    # Validate input
    if u_A isa AbstractVector{<:Real}
        @assert length(u_A) == 3 "u_A must be a 3-element vector"
    elseif u_A isa Function
        @assert applicable(u_A, 0) "function for u_A must be callable with one argument (t)"
        @assert length(u_A(0)) == 3 "u_A must be a function that returns a 3-element vector"
    end

    # Unpack unmodified kinematics of basis A
    @unpack p_A0,v_A,ω_A,vdot_A,ωdot_A = model.kinematicsBasisA

    # Set flags
    uAIsInput = true
    vAIsInput = false
    ωAIsInput = true
    vdotAIsInput = false
    ωdotAIsInput = false

    # Get new kinematics of basis A
    kinematicsBasisA = get_kinematics_basis_A(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A,uAIsInput,vAIsInput,ωAIsInput,vdotAIsInput,ωdotAIsInput)

    # Update initial states
    initialState = update_initial_states!(model.beams,model.assemblyData.elements,model.assemblyData.specialNodes,kinematicsBasisA,model.atmosphere)

    # Set on model
    setfield!(model, :kinematicsBasisA, kinematicsBasisA)
    setfield!(model, :initialState, initialState)

end
export update_u_A!


"""
    update_v_A!(model; v_A)

Updates the velocity of basis A (and related kinematics) in the model.

# Arguments
- `model::Model`: reference model

# Keyword Arguments
- `v_A::Union{AbstractVector{<:Real},<:Function}`: velocity of basis A
"""
function update_v_A!(model::Model; v_A::Union{AbstractVector{<:Real},<:Function})

    # Validate input
    if v_A isa AbstractVector{<:Real}
        @assert length(v_A) == 3 "v_A must be a 3-element vector"
    elseif v_A isa Function
        @assert applicable(v_A, 0) "function for v_A must be callable with one argument (t)"
        @assert length(v_A(0)) == 3 "v_A must be a function that returns a 3-element vector"
    end

    # Unpack unmodified kinematics of basis A
    @unpack p_A0,u_A,ω_A,vdot_A,ωdot_A = model.kinematicsBasisA

    # Set flags
    uAIsInput = false
    vAIsInput = true
    ωAIsInput = true
    vdotAIsInput = false
    ωdotAIsInput = false

    # Get new kinematics of basis A
    kinematicsBasisA = get_kinematics_basis_A(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A,uAIsInput,vAIsInput,ωAIsInput,vdotAIsInput,ωdotAIsInput)

    # Update initial states
    initialState = update_initial_states!(model.beams,model.assemblyData.elements,model.assemblyData.specialNodes,kinematicsBasisA,model.atmosphere)

    # Set on model
    setfield!(model, :kinematicsBasisA, kinematicsBasisA)
    setfield!(model, :initialState, initialState)

end
export update_v_A!


"""
    update_vdot_A!(model; vdot_A)

Updates the acceleration of basis A (and related kinematics) in the model.

# Arguments
- `model::Model`: reference model

# Keyword Arguments
- `vdot_A::Union{AbstractVector{<:Real},<:Function}`: acceleration of basis A
"""
function update_vdot_A!(model::Model; vdot_A::Union{AbstractVector{<:Real},<:Function})

    # Validate input
    if vdot_A isa AbstractVector{<:Real}
        @assert length(vdot_A) == 3 "vdot_A must be a 3-element vector"
    elseif vdot_A isa Function
        @assert applicable(vdot_A, 0) "function for vdot_A must be callable with one argument (t)"
        @assert length(vdot_A(0)) == 3 "vdot_A must be a function that returns a 3-element vector"
    end

    # Unpack unmodified kinematics of basis A
    @unpack p_A0,u_A,v_A,ω_A,ωdot_A = model.kinematicsBasisA

    # Set flags
    uAIsInput = false
    vAIsInput = false
    ωAIsInput = true
    vdotAIsInput = true
    ωdotAIsInput = false

    # Get new kinematics of basis A
    kinematicsBasisA = get_kinematics_basis_A(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A,uAIsInput,vAIsInput,ωAIsInput,vdotAIsInput,ωdotAIsInput)

    # Update initial states
    initialState = update_initial_states!(model.beams,model.assemblyData.elements,model.assemblyData.specialNodes,kinematicsBasisA,model.atmosphere)

    # Set on model
    setfield!(model, :kinematicsBasisA, kinematicsBasisA)
    setfield!(model, :initialState, initialState)

end
export update_vdot_A!


"""
    update_ω_A!(model; ω_A)

Updates the angular velocity of basis A (and related kinematics) in the model.

# Arguments
- `model::Model`: reference model

# Keyword Arguments
- `ω_A::Union{AbstractVector{<:Real},<:Function}`: angular velocity of basis A
"""
function update_ω_A!(model::Model; ω_A::Union{AbstractVector{<:Real},<:Function})

    # Validate input
    if ω_A isa AbstractVector{<:Real}
        @assert length(ω_A) == 3 "ω_A must be a 3-element vector"
    elseif ω_A isa Function
        @assert applicable(ω_A, 0) "function for ω_A must be callable with one argument (t)"
        @assert length(ω_A(0)) == 3 "ω_A must be a function that returns a 3-element vector"
    end

    # Unpack unmodified kinematics of basis A
    @unpack p_A0,u_A,v_A,vdot_A,ωdot_A = model.kinematicsBasisA

    # Set flags
    uAIsInput = true
    vAIsInput = false
    ωAIsInput = true
    vdotAIsInput = false
    ωdotAIsInput = false

    # Get new kinematics of basis A
    kinematicsBasisA = get_kinematics_basis_A(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A,uAIsInput,vAIsInput,ωAIsInput,vdotAIsInput,ωdotAIsInput)

    # Update initial states
    initialState = update_initial_states!(model.beams,model.assemblyData.elements,model.assemblyData.specialNodes,kinematicsBasisA,model.atmosphere)

    # Set on model
    setfield!(model, :kinematicsBasisA, kinematicsBasisA)
    setfield!(model, :initialState, initialState)

end
export update_ω_A!


"""
    update_ωdot_A!(model; ωdot_A)

Updates the angular acceleration of basis A (and related kinematics) in the model.

# Arguments
- `model::Model`: reference model

# Keyword Arguments
- `ωdot_A::Union{AbstractVector{<:Real},<:Function}`: angular acceleration of basis A
"""
function update_ωdot_A!(model::Model; ωdot_A::Union{AbstractVector{<:Real},<:Function})

    # Validate input
    if ωdot_A isa AbstractVector{<:Real}
        @assert length(ωdot_A) == 3 "ωdot_A must be a 3-element vector"
    elseif ωdot_A isa Function
        @assert applicable(ωdot_A, 0) "ωdot_A must be a function that can be called with a single argument (t)"
        @assert length(ωdot_A(0)) == 3 "ωdot_A must be a function that returns a 3-element vector"
    end

    # Unpack unmodified kinematics of basis A
    @unpack p_A0,u_A,v_A,ω_A,vdot_A = model.kinematicsBasisA

    # Set flags
    uAIsInput = true
    vAIsInput = false
    ωAIsInput = false
    vdotAIsInput = false
    ωdotAIsInput = true

    # Get new kinematics of basis A
    kinematicsBasisA = get_kinematics_basis_A(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A,uAIsInput,vAIsInput,ωAIsInput,vdotAIsInput,ωdotAIsInput)

    # Update initial states
    initialState = update_initial_states!(model.beams,model.assemblyData.elements,model.assemblyData.specialNodes,kinematicsBasisA,model.atmosphere)

    # Set on model
    setfield!(model, :kinematicsBasisA, kinematicsBasisA)
    setfield!(model, :initialState, initialState)

end
export update_ωdot_A!


"""
    add_BCs_to_model!(model; BCs)

Adds boundary conditions of the model.

# Arguments
- `model::Model`: reference model

# Keyword Arguments
- `BCs::Vector{BC}`: boundary conditions to add
"""
function add_BCs_to_model!(model::Model; BCs::Vector{BC})

    # Set new BCs on model
    push!(model.BCs, BCs...)

    # Update model assembly
    update_model_assembly!(model)

end
export add_BCs_to_model!


"""
    update_model_BCs!(model; BCs)

Updates the boundary conditions of the model.

# Arguments
- `model::Model`: reference model

# Keyword Arguments
- `BCs::Vector{BC}`: new array of boundary conditions
"""
function update_model_BCs!(model::Model; BCs::Vector{BC})

    # Set new BCs on model
    setfield!(model, :BCs, BCs)

    # Update model assembly
    update_model_assembly!(model)

end
export update_model_BCs!


"""
    update_model_links!(model; trimLoadsLinks,flapLinks)

Updates trim loads and flap deflection links of the model.

# Arguments
- `model::Model`: reference model

# Keyword Arguments
- `trimLoadsLinks::Vector{TrimLoadsLink}`: new array of trim loads links
- `flapLinks::Vector{FlapLink}`: new array of flap deflection links
"""
function update_model_links!(model::Model; trimLoadsLinks::Vector{TrimLoadsLink}=Vector{TrimLoadsLink}(),flapLinks::Vector{FlapLink}=Vector{FlapLink}())

    # Set new links on model
    setfield!(model, :trimLoadsLinks, trimLoadsLinks)
    setfield!(model, :flapLinks, flapLinks)

    # Update model assembly
    update_model_assembly!(model)

end
export update_model_links!


"""
    update_model_hinge_axis_constraints!(model; hingeAxisConstraints)

Updates the hinge axis constraint of the model.

# Arguments
- `model::Model`: reference model

# Keyword Arguments
- `hingeAxisConstraints::Vector{HingeAxisConstraint}`: new array of hinge axis constraints
"""
function update_model_hinge_axis_constraints!(model::Model; hingeAxisConstraints::Vector{HingeAxisConstraint})

    # Set new hinge axis constraints on model
    setfield!(model, :hingeAxisConstraints, hingeAxisConstraints)

    # Update model assembly
    update_model_assembly!(model)

end
export update_model_hinge_axis_constraints!


"""
    update_beam_orientation!(model; beam, p0)

Updates the assembly data of the model, after changing the initial orientation of a beam

# Arguments
- `model::Model`: reference model

# Keyword Arguments
- `beam::Beam`: the beam to update
- `p0::AbstractVector{<:Real}`: new initial Wiener-Milenkovic rotation parameters
"""
function update_beam_orientation!(model::Model; beam::Beam,p0::AbstractVector{<:Real})

    # Update beam
    update_beam_orientation!(beam=beam,p0=p0)

    # Update model
    update_model_assembly!(model)

end
export update_beam_orientation!


"""
    add_initial_conditions_to_beam!(model; beam,conditions)

Updates the initial states of a model's elements and nodes, after adding initial generalized displacements and velocities to a beam of the model

# Arguments
- `model::Model`: reference model

# Keyword arguments
- `beam::Beam`: reference beam
- `conditions::Dict{Symbol,<:Function}`: dictionary whose keys are initial condition types and values are the respective functions over arclength
"""
function add_initial_conditions_to_beam!(model::Model; beam::Beam,conditions::Dict{Symbol,<:Function})

    # Add initial conditions to the beam
    add_initial_conditions_to_beam!(beam,conditions=conditions)

    # Update initial state
    initialState = update_initial_states!(model.beams,model.assemblyData.elements,model.assemblyData.specialNodes,model.kinematicsBasisA,model.atmosphere)
    setfield!(model, :initialState, initialState)

end
export add_initial_conditions_to_beam!


"""
    add_point_inertia_to_beam!(model; beam,elementID,pointInertia)

Updates the model after adding point inertias to the beam

# Arguments
- `model::Model`: reference model

# Keyword arguments
- `beam::Beam`: reference beam
- `elementID::Int`: local ID (on the beam) of the element to which the point inertia is attached
- `pointInertia::PointInertia`: point inertia
"""
function add_point_inertia_to_beam!(model::Model; beam::Beam,elementID::Int,pointInertia::PointInertia)

    # Add point inertias to the beam
    add_point_inertia_to_beam!(beam=beam, elementID=elementID, pointInertia=pointInertia)

    # Update model
    update_model_assembly!(model)

end
export add_point_inertia_to_beam!


"""
    add_spring_to_beam!(model; beam,spring,nodeID,Δu0,Δp0,p)

Updates the model after adding simply-attached springs to a beam

# Arguments
- `model::Model`: reference model

# Keyword arguments
- `beam::Beam`: reference beam
- `spring::Spring`: spring
- `nodeID::Int`: local ID (on the beam) of the node to which the spring is attached
- `Δu0::Real`: translational pre-deformation
- `Δp0::Real`: rotational pre-deformation
- `p::AbstractVector{<:Real}`: Euler 3-2-1 rotation parameters from basis b of the attachment node to the spring's basis s
"""
function add_spring_to_beam!(model::Model; beam::Beam,spring::Spring, nodeID::Int,Δu0::Real=0,Δp0::Real=0,p::AbstractVector{<:Real}=zeros(3))

    # Add spring to the beam
    add_spring_to_beam!(beam=beam, spring=spring, nodeID=nodeID, Δu0=Δu0, Δp0=Δp0, p=p)

    # Update model
    update_model_assembly!(model)

end
export add_spring_to_beam!


"""
    add_spring_to_beams!(model; beams, spring, nodesLocalIDs, Δu0, Δp0, p)

Updates the model after adding a doubly-attached spring to the beams

# Arguments
- `model::Model`: reference model

# Keyword arguments
- `beams::Vector{Beam}`: reference beams
- `spring::Spring`: doubly-attached spring
- `nodesLocalIDs`: local IDs (on the respective beams) of the nodes to which the spring is attached
- `Δu0::Real`: translational pre-deformation
- `Δp0::Real`: rotational pre-deformation
- `p::AbstractVector{<:Real}`: Euler 3-2-1 rotation parameters from basis b of the attachment node to the spring's basis s
"""
function add_spring_to_beams!(model::Model; beams::Vector{Beam},spring::Spring,nodesLocalIDs::AbstractVector{Int},Δu0::Real=0,Δp0::Real=0,p::AbstractVector{<:Real}=zeros(3))

    # Add spring to the beams
    add_spring_to_beams!(beams=beams, spring=spring, nodesLocalIDs=nodesLocalIDs, Δu0=Δu0, Δp0=Δp0, p=p)

    # Update model
    update_model_assembly!(model)

end
export add_spring_to_beams!


"""
    add_surface_to_beam!(model; beam,surface)

Updates the model after adding an aerodynamic surface to the beam

# Arguments
- `model::Model`: reference model

# Keyword arguments
- `beam::Beam`: reference beam
- `surface::AeroSurface`: aerodynamic surface
"""
function add_surface_to_beam!(model::Model; beam::Beam,surface::AeroSurface)

    # Add aerodynamic surface to the beam
    add_surface_to_beam!(beam,surface=surface)  

    # Update model
    update_model_assembly!(model)

end
export add_surface_to_beam!


"""
    remove_all_initial_conditions_from_beams!(model; beams)

Updates the model after removing all initial conditions from the beams

# Arguments
- `model::Model`: reference model

# Keyword arguments
- `beams::AbstractVector{Beam}`: array of beams
"""
function remove_all_initial_conditions_from_beams!(model::Model; beams::AbstractVector{Beam})

    # Remove all initial conditions from the beams
    remove_all_initial_conditions_from_beams!(beams=beams)

    # Update initial state
    initialState = update_initial_states!(model.beams,model.assemblyData.elements,model.assemblyData.specialNodes,model.kinematicsBasisA,model.atmosphere)
    setfield!(model, :initialState, initialState)

end
export remove_all_initial_conditions_from_beams!


"""
    remove_all_point_inertias_from_beams!(model; beams=beams)

Updates the model after removing all point inertias attached to the beams

# Arguments
- `model::Model`: reference model

# Keyword arguments
- `beams::AbstractVector{Beam}`: array of beams
"""
function remove_all_point_inertias_from_beams!(model::Model; beams::AbstractVector{Beam})

    # Remove all point inertias from the beams
    remove_all_point_inertias_from_beams!(beams=beams)

    # Update model
    update_model_assembly!(model)

end
export remove_all_point_inertias_from_beams!


"""
    remove_all_springs_from_beams!(model; beams)

Updates the model after removing all springs attached to the beams

# Arguments
- `model::Model`: reference model

# Keyword arguments
- `beams::AbstractVector{Beam}`: array of beams
"""
function remove_all_springs_from_beams!(model::Model; beams::AbstractVector{Beam})

    # Remove all springs from the beams
    remove_all_springs_from_beams!(beams=beams)

    # Update model
    update_model_assembly!(model)

end
export remove_all_springs_from_beams!



"""
    remove_surface_from_beams!(model; beams)

Updates the model after removing the aerodynamic surface attached to the beams

# Arguments
- `model::Model`: reference model

# Keyword arguments
- `beams::AbstractVector{Beam}`: array of beams
"""
function remove_surface_from_beams!(model::Model; beams::AbstractVector{Beam})

    # Remove aerodynamic surface from the beams
    remove_surface_from_beams!(beams=beams)

    # Update model
    update_model_assembly!(model)

end
export remove_surface_from_beams!


"""
    remove_gust!(model)
    
Removes the gust from the model and updates the model assembly

# Arguments
- `model::Model`: reference model
"""
function remove_gust!(model::Model)

    # Remove gust from the model
    setfield!(model, :gust, nothing)

    # Update model assembly
    update_model_assembly!(model)

end
export remove_gust!


## PRIVATE API
# ------------------------------------------------------------------------------
# Validates the inputs to the model and sets correct types
function validate_model(beams,firstNodePosition,gravityVector,connections,p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A)

    # Validate beams
    @assert !isempty(beams) "model cannot be empty of beams"

    # Validate initial position of the first node of the model
    @assert length(firstNodePosition) == 3 "firstNodePosition must have length 3"
    firstNodePosition = SVector{3,Float64}(firstNodePosition)

    # Validate gravity vector
    @assert length(gravityVector) == 3 "gravityVector must have length 3"
    gravityVector = SVector{3,Float64}(gravityVector)

    # Validate connections
    validate_connections(beams,connections)

    # Valide motion of basis A
    p_A0,uAIsInput,vAIsInput,ωAIsInput,vdotAIsInput,ωdotAIsInput = validate_motion_of_basis_A(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A)

    return firstNodePosition,gravityVector,connections,p_A0,uAIsInput,vAIsInput,ωAIsInput,vdotAIsInput,ωdotAIsInput

end


# Validates beam connections
function validate_connections(beams,connections)

    # Initialize a dictionary of beam IDs for validation
    beamIDs = Dict{Beam,Nothing}(beam => nothing for beam in beams)

    # Loop connections and validate each
    for connection in connections
        @assert connection isa Tuple && length(connection) == 2 "each connection must be a tuple of the form (beamConnection, nodeConnection)"
        beamConnection = connection[1]
        nodeConnection = connection[2]
        @assert beamConnection isa Pair{Beam,Beam} "each beamConnection must be a Pair{Beam,Beam}"
        @assert nodeConnection isa Pair{Int,Int} "each nodeConnection must be a Pair{Int,Int}"
        @assert haskey(beamIDs, beamConnection.first) && haskey(beamIDs, beamConnection.second) "all beams in connections must be included in the model beams array"
    end

    return connections

end

# Validates the motion of basis A, sets correct types, and gets flags for which motion variables are inputs
function validate_motion_of_basis_A(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A)

    # Validate initial rotation of basis A
    @assert length(p_A0) == 3 "p_A0 must have length 3"
    p_A0 = SVector{3,Float64}(p_A0)

    # Validate motion functions
    for x in [u_A,v_A,ω_A,vdot_A,ωdot_A]
        if !isnothing(x)
            @assert (x isa Function && applicable(x, 0) && length(x(0)) == 3) || (x isa AbstractVector{<:Real} && length(x) == 3) "motion variables must be either a 3-element vector or a function of one argument (t) that returns a 3-element vector"
        end
    end

    # Set flags for motion of basis A being inputs
    uAIsInput    = !isnothing(u_A)
    vAIsInput    = !isnothing(v_A)
    ωAIsInput    = !isnothing(ω_A)
    vdotAIsInput = !isnothing(vdot_A)
    ωdotAIsInput = !isnothing(ωdot_A)

    # Check that only one of u_A, v_A, vdot_A is specified and that only one of ω_A, ωdot_A is specified
    @assert count([uAIsInput, vAIsInput, vdotAIsInput]) ≤ 1 "Specify only one of u_A, v_A, vdot_A"
    @assert count([ωAIsInput, ωdotAIsInput]) ≤ 1 "Specify only one of ω_A, ωdot_A"

    return p_A0,uAIsInput,vAIsInput,ωAIsInput,vdotAIsInput,ωdotAIsInput

end


# Computes the motion variables of basis A
function get_kinematics_basis_A(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A,uAIsInput,vAIsInput,ωAIsInput,vdotAIsInput,ωdotAIsInput)

    # Set as functions of time
    u_A    = set_as_time_function(u_A)
    v_A    = set_as_time_function(v_A)
    ω_A    = set_as_time_function(ω_A)
    vdot_A = set_as_time_function(vdot_A)
    ωdot_A = set_as_time_function(ωdot_A)

    # Compute missing motion variables from given ones as functions of time
    if uAIsInput
        v_A    = t -> ForwardDiff.derivative(u_A, t)
        vdot_A = t -> ForwardDiff.derivative(v_A, t)
    elseif vAIsInput
        u_A    = t -> QuadGK.quadgk(v_A, 0, t)[1]
        vdot_A = t -> ForwardDiff.derivative(v_A, t)
    elseif vdotAIsInput
        v_A = t -> QuadGK.quadgk(vdot_A, 0, t)[1]
        u_A = t -> QuadGK.quadgk(v_A, 0, t)[1]
    end
    if ωAIsInput
        ωdot_A = t -> ForwardDiff.derivative(ω_A, t)
    elseif ωdotAIsInput
        ω_A = t -> QuadGK.quadgk(ωdot_A, 0, t)[1]
    end

    # Set kinematics object for basis A
    kinematicsBasisA = KinematicsBasisA(p_A0,u_A,v_A,ω_A,vdot_A,ωdot_A)

    return kinematicsBasisA

end


# Sets the motion variables of basis A as functions of time
function set_as_time_function(x)
    if isnothing(x)
        return t -> z3
    elseif x isa Function
        return t -> SVector{3}(x(t))
    else #if x isa AbstractVector{<:Real}
        constVal = SVector{3,Float64}(x)
        return t -> constVal
    end
end


# Assembles the beams into a model and gets related information for the model, including the assembly data of the beams, elements, special nodes, BCs, trim loads and flap links, and hinge axis constraints
function get_model_assembly_data!(beams,connections,firstNodePosition,trimLoadsLinks,flapLinks,hingeAxisConstraints,BCs,gust)

    # Update flap links assembly data
    update_flap_links_assembly_data!(flapLinks)

    # Set correct dimensions for aerodynamic vectors of elements
    redimensionalize_elemental_aero_data!(beams,gust)

    # Split connections into beam and node connections
    beamConnections,nodeConnections = split_connections(connections)

    # Assemble the beams' elements into a model
    elements,nElementsTotal,nNodesTotal,elementNodes,r_n = assemble_elements!(beams,beamConnections,nodeConnections,firstNodePosition)

    # Update the assembly data of the trim loads links
    update_load_links_assembly_data!(trimLoadsLinks)

    # Compute force scaling
    forceScaling = get_force_scaling(elements)

    # Compute initial inertia properties of the assembly
    initialInertiaProperties = get_initial_inertia_properties(elements)

    # Get BC'ed nodes of the model
    BCedNodes = get_BCed_nodes(BCs)

    # Create and get special nodes
    specialNodes,specialNodesGlobalIDs = get_special_nodes(beams,elements,nNodesTotal,elementNodes,BCs,BCedNodes)

    # Set system indices
    systemOrder,rangeNominal,rangeAeroStates = set_system_indices!(elements,nNodesTotal,specialNodes,specialNodesGlobalIDs,hingeAxisConstraints)

    # Update double attachment springs data for all special nodes
    update_double_attachment_springs_data!(specialNodes)

    # Update trim assembly data and get number of trim variables (trim loads and flap deflections)
    nTrimVariables,nTrimLoads,nTrimδ,rangeTrimLoads,rangeTrimδ = update_trim_assembly_data!(beams,specialNodes,trimLoadsLinks,flapLinks,systemOrder)

    # Create assembly data
    assemblyData = ModelAssemblyData(
        elements,
        group_elements_by_type(elements),
        nElementsTotal,
        nNodesTotal,
        elementNodes,
        r_n,
        forceScaling,
        1/forceScaling,
        BCedNodes,
        specialNodes,
        specialNodesGlobalIDs,
        systemOrder,
        nTrimVariables,
        nTrimLoads,
        nTrimδ,
        rangeNominal,
        rangeAeroStates,
        rangeTrimLoads,
        rangeTrimδ,
        initialInertiaProperties
    )

    return assemblyData,BCs

end


# Groups the elements into a tuple of concretely-typed vectors, which share the elements with the original vector. The groups are runs of consecutive elements of the same concrete type, so that looping over the groups visits the elements in the very same (assembly) order, which matters: the equations' indices and the "equations already set" flags of the elements' assembly indices are defined by that order
function group_elements_by_type(elements)

    # Maximum number of groups (the tuple type grows with it)
    maxGroups = 16

    # Initial and final indices of the runs of elements with the same type
    types = typeof.(elements)
    firsts = [1; findall(i -> types[i] != types[i-1], 2:length(elements)) .+ 1]
    lasts = [firsts[2:end] .- 1; length(elements)]

    # Fall back to a single (abstractly-typed) group if there are too many runs
    length(firsts) > maxGroups && return (elements,)

    return Tuple(convert(Vector{types[a]}, elements[a:b]) for (a,b) in zip(firsts,lasts))

end


# Updates the assembly data of the flap links
function update_flap_links_assembly_data!(flapLinks)

    # Loop flap links
    for flapLink in flapLinks
        @unpack masterBeam,slaveBeams,δMultipliers = flapLink
        @unpack δFun,δdotFun,δddotFun,δIsTrimVariable = masterBeam.attachments.surface.assemblyData
        # Loop slave beams
        slaveSurfaces = [slaveBeam.attachments.surface for slaveBeam in slaveBeams]
        for (i,slaveSurface) in enumerate(slaveSurfaces)
            # Update slave surface
            setfield!(slaveSurface, :assemblyData, AeroSurfaceAssemblyData(slaveSurface, δFun, δdotFun, δddotFun, δMultipliers[i], δIsTrimVariable))
        end
    end

end


# Sets the correct dimensions for aerodynamic vectors in elemental fields, now that the gust and flap trim states are known
function redimensionalize_elemental_aero_data!(beams,gust)

    # Set flag
    withGust = !isnothing(gust)

    # Loop over beams
    for beam in beams
        # Unpack beam data
        @unpack surface = beam.attachments
        @unpack elements = beam
        # Skip if no aerodynamic surface is attached to the beam
        isnothing(surface) && continue
        # Unpack surface data
        @unpack nPitchPlungeStates,nFlapStates,nGustStates = surface.stateRanges
        @unpack δIsTrimVariable = surface.assemblyData
        # Loop over elements and set correct dimensions for aerodynamic vectors
        for (e,element) in enumerate(elements)
            # Unpack element aero data
            @unpack isFlapped = element.aero.geometry
            # Element's number of aerodynamic states
            N = nPitchPlungeStates + nFlapStates * isFlapped + nGustStates * !isnothing(gust)
            # Element's number of trim flap states
            Nδ = isFlapped && δIsTrimVariable ? 
                 1 : 
                 0
            # Rebuild element with correct dimensions for aerodynamic arrays
            elements[e] = rebuild_element(element,N,Nδ,withGust)
        end
    end

end


# Splits combined connections into beam and node vectors
function split_connections(connections)

    beamConnections = Pair{Beam,Beam}[connection for (connection,_) in connections]
    nodeConnections = Pair{Int,Int}[nodes for (_,nodes) in connections]

    return beamConnections,nodeConnections

end


# Assembles the elements of the beams, gets related information for the model, and updates beams' (and their elements') assembly data
function assemble_elements!(beams,beamConnections,nodeConnections,firstNodePosition)

    # Initialize assembly arrays
    elements = Vector{Element}()
    elementNodes = Vector{Vector{Int}}()
    r_n = Vector{Vector{Float64}}()

    # Initialize element and node counts in the assembly
    elementCount = 0
    nodeCount = 0

    # Common rounding tolerance for the assembly coordinates (must be the same for all beams, so that coincident nodes of different beams are rounded identically)
    roundTol = 1e-8*minimum(element.geometry.Δℓ for beam in beams for element in beam.elements)

    # Loop over beams
    for (ID,beam) in enumerate(beams)

        @unpack nElements = beam.discretization

        # Add current beam's elements to model's array
        for element in beam.elements
            push!(elements,element)
        end

        # Set range of global IDs of the elements of the current beam
        beamElementsGlobalIDs = collect(1:nElements) .+ elementCount

        # Increment element count
        elementCount += nElements

        # Get beam connections
        connectedBeamsIDs,connectedNodesThisLocalID,connectedNodesOtherLocalID = get_current_beam_connections_IDs(beams,ID,beamConnections,nodeConnections)

        # Update elements' nodes global IDs assembly array (all elements in the model)
        elementNodes = update_element_nodes(elementNodes,beams,ID,nElements,connectedBeamsIDs,connectedNodesThisLocalID,connectedNodesOtherLocalID,nodeCount,beamElementsGlobalIDs)

        # Array of global nodal IDs for each of the current beam's elements
        elementNodesGlobalIDs = [elementNodes[e] for e in beamElementsGlobalIDs]

        # Array of global nodal IDs for the current beam's nodes
        beamNodesGlobalIDs = unique(vcat(elementNodes[beamElementsGlobalIDs]...))

        # Update node count
        nodeCount = maximum(beamNodesGlobalIDs)

        # Update assembly coordinates
        r_n,r_nBeam,rElem,r_n1Elem,r_n2Elem = update_assembly_coordinates(r_n,beam,ID,nElements,beamNodesGlobalIDs,firstNodePosition,roundTol)

        # Update beam assembly data
        setfield!(beam, :assemblyData, BeamAssemblyData(ID,beamElementsGlobalIDs,beamNodesGlobalIDs,r_nBeam,connectedBeamsIDs,connectedNodesThisLocalID,connectedNodesOtherLocalID))

        # Update beam's elements' assembly data
        for (e, element) in enumerate(beam.elements)
            setfield!(element, :assemblyData, ElementAssemblyData(element; parentBeamID=ID, globalID=beamElementsGlobalIDs[e], nodesGlobalIDs=elementNodesGlobalIDs[e], r=rElem[e], r_n1=r_n1Elem[e], r_n2=r_n2Elem[e]))
        end

    end

    # Total number of elements in the assembly
    nElementsTotal = elementCount

    # Total number of nodes in the assembly
    nNodesTotal = nodeCount

    # Set correct types
    elementNodes = [SVector{2,Int}(x) for x in elementNodes]
    r_n = [SVector{3,Float64}(x) for x in r_n]

    return elements,nElementsTotal,nNodesTotal,elementNodes,r_n

end


# Gets the connected beams and the nodes of the current and connected beams that are connected
function get_current_beam_connections_IDs(beams,ID,beamConnections,nodeConnections)

    # Initialize arrays
    connectedBeamsIDs = Vector{Int}()
    connectedNodesThisLocalID = Vector{Int}()
    connectedNodesOtherLocalID = Vector{Int}()

    # Current beam
    thisBeam = beams[ID]

    # Loop over beam connections
    for (connection,nodes) in zip(beamConnections,nodeConnections)
        # Unpack beams and nodes of the connection
        beam1 = connection.first
        beam2 = connection.second
        node1 = nodes.first
        node2 = nodes.second
        # Set connected beam ID and nodes of the current and connected beams
        if beam1 === thisBeam
            connectedBeamID = findfirst(==(beam2), beams)
            nodeThis = node1
            nodeOther = node2
        else
            continue
        end
        push!(connectedBeamsIDs, connectedBeamID)
        push!(connectedNodesThisLocalID, nodeThis)
        push!(connectedNodesOtherLocalID, nodeOther)
    end

    return connectedBeamsIDs,connectedNodesThisLocalID,connectedNodesOtherLocalID

end


# Updates the array of global IDs of the nodes of the elements of a beam in the assembly
function update_element_nodes(elementNodes,beams,ID,nElements,connectedBeamsIDs,connectedNodesThisLocalID,connectedNodesOtherLocalID,nodeCount,beamElementsGlobalIDs)

    # Increment elements' nodes' global ID array assuming current beam is simply connected
    for e in beamElementsGlobalIDs
        push!(elementNodes, [e, e+1])
    end

    # Update elements' nodes' global ID array in case the beam is not simply connected
    if ID > 1 && elementNodes[beamElementsGlobalIDs][1][2] > nodeCount+1
        diff = elementNodes[beamElementsGlobalIDs][1][2] - (nodeCount+1)
        elementNodes[beamElementsGlobalIDs] .= [x .- [diff,diff] for x in elementNodes[beamElementsGlobalIDs]]
    end
    for (connectedBeamID,connectedNodeThisID,connectedNodeOtherID) in zip( connectedBeamsIDs,connectedNodesThisLocalID,connectedNodesOtherLocalID)
        # Connected beam
        connectedBeam = beams[connectedBeamID]
        # Local IDs of the elements of the current beam to which the connected node belongs
        currentBeamConnectedElementsLocalID = connectedNodeThisID == 1 ? 
            [1] : 
            connectedNodeThisID == nElements+1 ? 
                [nElements] : 
                [connectedNodeThisID-1; connectedNodeThisID]
        # Global IDs of the above
        currentBeamConnectedElementsGlobalID = beamElementsGlobalIDs[currentBeamConnectedElementsLocalID]
        # Sides of the connected node on the current beam's connected elements
        currentBeamNodeSides = connectedNodeThisID == 1 ? 
            [1] : 
            connectedNodeThisID == nElements+1 ? 
                [2] : 
                [2; 1]
        # Global ID of the node on the connected beam
        connectedBeamsConnectedNodeGlobalID = connectedBeam.assemblyData.nodesGlobalIDs[connectedNodeOtherID]
        # In case there is more than one connected element in the current beam, increase the indices of the nodes up to the connected node
        if length(currentBeamConnectedElementsGlobalID) > 1
            elementNodes[beamElementsGlobalIDs[1:currentBeamConnectedElementsLocalID[1]]] .= [x .+ [1,1] for x in elementNodes[beamElementsGlobalIDs[1:currentBeamConnectedElementsLocalID[1]]]]
        end
        # Loop current beam's connected elements
        for (currentBeamConnectedElementGlobalID,currentBeamNodeSide) in zip(currentBeamConnectedElementsGlobalID,currentBeamNodeSides)
            # Update connected element's nodes' global ID
            elementNodes[currentBeamConnectedElementGlobalID][currentBeamNodeSide] = connectedBeamsConnectedNodeGlobalID
        end
    end

    return elementNodes
    
end


# Sets the coordinates of the current beam and its elements and updates the assembly coordinates
function update_assembly_coordinates(r_n,beam,ID,nElements,beamNodesGlobalIDs,firstNodePosition,roundTol)

    # Get the global ID of the first node of the current beam
    firstNodeGlobalID = minimum(beamNodesGlobalIDs[1:end-1])

    # Initialize array of nodal coordinates of the current beam's elements
    r_nBeam = Vector{Vector{Float64}}(undef, nElements+1)

    # Initialize elements' nodal and midpoint coordinates arrays
    rElem = Vector{Vector{Float64}}(undef,nElements)
    r_n1Elem = Vector{Vector{Float64}}(undef,nElements)
    r_n2Elem = Vector{Vector{Float64}}(undef,nElements)
    
    # Set coordinates (defined in basis A) at beam and element levels
    for (e,element) in enumerate(beam.elements)
        @unpack k,R0,firstNodeOffset = beam.geometry
        @unpack x1,x1_n1,x1_n2 = element.geometry
        # Initialize element's coordinates as those relative to beam origin
        r = position_vector(R0,k,x1)
        r_n1 = position_vector(R0,k,x1_n1)
        r_n2 = position_vector(R0,k,x1_n2)
        # Increment with beam's first node position offset
        r .+= firstNodeOffset
        r_n1 .+= firstNodeOffset
        r_n2 .+= firstNodeOffset
        # If in the first beam, add the position of the model's first node to element's coordinates, otherwise, add coordinates of the beam's first node
        rFirstNode = ID == 1 ? 
                     firstNodePosition : 
                     r_n[firstNodeGlobalID]
        r .+= rFirstNode
        r_n1 .+= rFirstNode
        r_n2 .+= rFirstNode
        # Deal with floating point errors
        round_off!(r, roundTol)
        round_off!(r_n1, roundTol)
        round_off!(r_n2, roundTol)
        # Deal with signed zeros
        r = map(x -> iszero(x) ? 0.0 : x, r)
        r_n1 = map(x -> iszero(x) ? 0.0 : x, r_n1)
        r_n2 = map(x -> iszero(x) ? 0.0 : x, r_n2)
        # Set current element's coordinates in beam and element arrays
        rElem[e] = r
        r_n1Elem[e] = r_n1
        r_n2Elem[e] = r_n2
        if e == 1
            r_nBeam[e] = r_n1
            r_nBeam[e+1] = r_n2
        else
            r_nBeam[e+1] = r_n2
        end
    end

    # Update assembly coordinates (add current beam's coordinates and remove duplicates from repeated nodes)
    push!(r_n,r_nBeam...)
    unique!(r_n)

    # Set correct types for coordinates
    r_nBeam = [SVector{3}(x) for x in r_nBeam]
    rElem = [SVector{3}(x) for x in rElem]
    r_n1Elem = [SVector{3}(x) for x in r_n1Elem]
    r_n2Elem = [SVector{3}(x) for x in r_n2Elem]

    return r_n,r_nBeam,rElem,r_n1Elem,r_n2Elem

end


# Updates the assembly data of the trim loads links
function update_load_links_assembly_data!(trimLoadsLinks)

    # Loop trim load links
    for loadsLink in trimLoadsLinks
        # Unpack trim loads link data
        @unpack masterBC,slaveBCs,masterBeam,slaveBeams = loadsLink.masterAndSlave
        # Local IDs of the master and slave nodes
        masterNodeLocalID = masterBC.definition.node
        slaveNodesLocalIDs = [BC.definition.node for BC in slaveBCs]
        # Master node global IDs
        masterNodeGlobalID = masterBeam.assemblyData.nodesGlobalIDs[masterNodeLocalID]
        # Slave node global IDs
        slaveNodesGlobalIDs = Vector{Int}()
        for (slaveBeam,slaveNodeLocalID) in zip(slaveBeams,slaveNodesLocalIDs)
            slaveNodeGlobalID = slaveBeam.assemblyData.nodesGlobalIDs[slaveNodeLocalID]
            push!(slaveNodesGlobalIDs,slaveNodeGlobalID)
        end
        # Set assembly data of the trim loads link
        setfield!(loadsLink, :assemblyData, TrimLoadsLinkAssemblyData(masterNodeGlobalID,slaveNodesGlobalIDs))
    end

end


# Computes an appropriate force scaling for the linear system of equations
function get_force_scaling(elements)

    # Initialize
    nNonzero = 0
    sumNonzero = 0

    # Count nonzero compliance matrix entries and their sum
    for element in elements
        C = element.sectional.C
        for c in C
            if !iszero(c)
                nNonzero += 1
                sumNonzero += c
            end
        end
    end

    # Power factor (inverse of the average nonzero compliance matrix entry)
    α = nNonzero / sumNonzero

    # Force scaling factor
    forceScaling = nNonzero > 0 ? 
                  (α > 2 ? 
                    nextpow(2, α) : 
                    2 / nextpow(2, 1/α)
                  ) : 
                  1.0       

    return forceScaling

end


# Computes the inertia properties of the undeformed assembly
function get_initial_inertia_properties(elements)

    # Elements' mass array
    elementsMass = [e.geometry.Δℓ * e.sectional.μ for e in elements]

    # Elements position vector from the assembly's origin to the element's center of mass, resolved in basis A
    elementsR = [e.assemblyData.r + e.geometry.R0 * axial(e.sectional.ηtilde) for e in elements]

    # Elements' inertia matrices, resolved in basis A
    elementsI = [e.geometry.Δℓ * e.geometry.R0 * e.sectional.I[4:6,4:6] * e.geometry.R0' for e in elements]

    # Mass of the assembly
    mass = sum(elementsMass)
    
    # Position of the center of mass, resolved in basis A
    centerOfMass = zeros(3)
    for e in eachindex(elements)
        centerOfMass .+= elementsR[e] * elementsMass[e] / mass
    end

    # Inertia matrix of the assembly about its center of mass, resolved in basis A
    inertiaMatrix = zeros(3,3)
    for e in eachindex(elements)
        r = elementsR[e] - centerOfMass # position vector from the assembly's center of mass to element's center of mass, resolved in basis A
        inertiaMatrix .+= elementsI[e] + elementsMass[e] * (dot(r,r) * I3 - r*r')
    end

    # Convert to SVector and SMatrix types
    centerOfMass = SVector{3,Float64}(centerOfMass)
    inertiaMatrix = SMatrix{3,3,Float64,9}(inertiaMatrix)

    # Set initial inertia properties object
    initialInertiaProperties = InitialInertiaProperties(mass,centerOfMass,inertiaMatrix)

    return initialInertiaProperties

end


# Gets the BC'ed nodes
function get_BCed_nodes(BCs)

    # Initialize list of BC'ed nodes
    BCedNodes = Int[]

    # Loop BCs and add that BC's node to the BC'ed nodes list
    for BC in BCs
        push!(BCedNodes,BC.definition.beam.assemblyData.nodesGlobalIDs[BC.definition.node])
    end

    # Remove duplicates
    unique!(BCedNodes)

    return BCedNodes
end


# Gets the special nodes in the system of equations: connection, boundary, and BC'ed nodes
function get_special_nodes(beams,elements,nNodesTotal,elementNodes,BCs,BCedNodes)

    # Get the nodes that have springs attached to them and the list of springs on each node
    springedNodes,springsOnNodes = get_spring_assembly_data(beams,nNodesTotal)

    # Initialize array of special nodes
    specialNodes = Vector{SpecialNode}()

    # List of all nodes
    allNodesList = 1:nNodesTotal

    # Initialize flag for nodes being special
    isSpecial = falses(nNodesTotal)

    # BC'ed nodes (any essential or non-zero natural BC) are special 
    isSpecial[BCedNodes] .= true

    # Nodes with attached springs are special
    isSpecial[springedNodes] .= true

    # Initialize special nodes counter
    s = 0

    # Loop over nodes
    for node in allNodesList
        # Elements connected to this node and their IDs
        connectedElementsGlobalIDs = findall(x -> node in x, elementNodes)
        connectedElements = elements[connectedElementsGlobalIDs]
        # Number of elements that are connected to the node
        nConnectedElements = length(connectedElements)
        # If this is either a boundary (connected to only 1 element) or a connection node (connected to 3 or more elements), set as special 
        if nConnectedElements != 2                         
            isSpecial[node] = true          
        end
        # Skip non-special nodes
        if !isSpecial[node]
            continue
        end
        # Increment special nodes count
        s += 1
        # Find in which side and ζ position of the elements that node is 
        sideOnElements = Vector{Int}(undef, nConnectedElements)
        ζonElements = Vector{Int}(undef, nConnectedElements)
        for (i, e) in enumerate(connectedElementsGlobalIDs)
            sideOnElements[i] = elementNodes[e][1] == node ? 
                                1 :
                                2
            ζonElements[i] = elementNodes[e][1] == node ?
                             -1 :
                             1
        end
        # Get the global ID of the node
        globalID = connectedElements[1].assemblyData.nodesGlobalIDs[sideOnElements[1]]
        # Initialize array of BCs on that node
        nodesBCs = Vector{BC}()
        # Add BCs for BC'ed nodes
        for BC in BCs
            if node == BC.definition.beam.assemblyData.nodesGlobalIDs[BC.definition.node]
                push!(nodesBCs,BC)
            end
        end
        # Add special node to array
        push!(specialNodes, create_SpecialNode(globalID,connectedElementsGlobalIDs,ζonElements,nodesBCs,springsOnNodes[node]))
    end

    # Global IDs of the special nodes
    specialNodesGlobalIDs = allNodesList[isSpecial]

    # Global node to special node ID dictionary
    global2special = Dict(g => i for (i,g) in enumerate(specialNodesGlobalIDs))

    # Loop over special nodes
    for specialNode in specialNodes
        # Loop over node's springs
        for spring in specialNode.springs
            # Unpack spring data
            @unpack nodesGlobalIDs = spring.assemblyData
            # Get the special IDs of the spring's nodes
            nodesSpecialIDs = [global2special[nodesGlobalIDs[i]] for i in eachindex(nodesGlobalIDs)]
            # Update
            setfield!(spring, :assemblyData, SpringAssemblyData(spring,nodesGlobalIDs,nodesSpecialIDs))
        end
    end

    return specialNodes,specialNodesGlobalIDs

end


# Gets the list of springed nodes and the list of springs on each node
function get_spring_assembly_data(beams,nNodesTotal)

    # Initialize list of springed nodes and list of springs on nodes
    springedNodes = Int[]
    springsOnNodes = [Spring[] for _ in 1:nNodesTotal]

    # Collect all springs from beams
    springs = unique([spring for beam in beams for spring in beam.attachments.springs])

    # Loop springs
    for spring in springs
        # Unpack
        @unpack beamConstructionIDs,elemsLocalIDs,nodesSides = spring.assemblyData
        # Initialize global IDs of the spring's nodes
        nodesGlobalIDs = Vector{Int}(undef,length(nodesSides))
        for (i,constructionID,e,side) in zip(1:length(beamConstructionIDs),beamConstructionIDs,elemsLocalIDs,nodesSides)
            # Find the beam with the given construction ID
            beam = beams[findfirst(b -> b.constructionID == constructionID, beams)]
            # Global ID of current attachment node
            ID = beam.elements[e].assemblyData.nodesGlobalIDs[side]
            # Update spring's node global ID
            nodesGlobalIDs[i] = ID
            # Update list of springed nodes and list of springs on nodes
            push!(springedNodes, ID)
            push!(springsOnNodes[ID], spring)
        end
        # Update spring's assembly data
        setfield!(spring, :assemblyData, SpringAssemblyData(spring,nodesGlobalIDs))
    end

    return springedNodes,springsOnNodes
    
end


# Gets the indices (for equations and DOFs) of the system of equations
function set_system_indices!(elements,nNodesTotal,specialNodes,specialNodesGlobalIDs,hingeAxisConstraints)
    
    ## Initialization
    #---------------------------------------------------------------------------
    # Flag indicating whether indices of a node's equations/states indices have been assigned
    assigned = falses(nNodesTotal)
    # Array with element on which node's equations have been assigned (first column contains the elements, and second column which local node it is in that element)
    e_assigned = Vector{Vector{Int}}()
    # Counters for equations/states
    i_equations = 1                   
    i_states = 1
    # List of all aerodynamic DOFs
    rangeAeroStates = Int[]
    
    ## Set indices of elemental and standard nodes' structural equations/states, and DOFs of special nodes
    #---------------------------------------------------------------------------
    # Loop over elements
    for (e, element) in enumerate(elements)

        # Unpack element data
        @unpack nodesGlobalIDs = element.assemblyData
        nAeroStates = length(element.state.midpoint.χ)

        # --- First node --- #
        # ----------------------------------------------------------------------

        # Get element's first node's global ID
        n = nodesGlobalIDs[1]

        # Initialize flags for first node
        isSpecialNode1 = false
        eqsNode1Set = false

        # Check if the node's equations/states have been assigned indices
        if !assigned[n]
            # Update assignment flag and pointer to element
            assigned[n] = true
            push!(e_assigned,[e, 1])
            # Indices for equilibrium and compatibility equations of element's first node
            eqs_Fu1 = SVector{3,Int}(i_equations+0:i_equations+2)
            eqs_Fp1 = SVector{3,Int}(i_equations+3:i_equations+5)
            eqs_FF1 = SVector{3,Int}(i_equations+6:i_equations+8)
            eqs_FM1 = SVector{3,Int}(i_equations+9:i_equations+11)
            # Update equations' counter
            i_equations += 12
            # Check whether this is a special node - has nodal states
            if n in specialNodesGlobalIDs
                # Get corresponding special node
                s = findfirst(x -> n in x, specialNodesGlobalIDs)
                specialNode = specialNodes[s]
                # Indices for nodal states
                DOF_uF = SVector{3,Int}(i_states+0:i_states+2)
                DOF_pM = SVector{3,Int}(i_states+3:i_states+5)
                # Update states' counter
                i_states += 6
                # Set flag
                isSpecialNode1 = true
                # Update special node's assembly data
                setfield!(specialNode, :assemblyData, SpecialNodeAssemblyData(specialNode, DOF_uF, DOF_pM))
            end
        else   
            # Update flag for equations having already been set
            eqsNode1Set = true
            # Set indices for equilibrium equations as those already assigned
            if e_assigned[n][2] == 1
                eqs_Fu1 = elements[e_assigned[n][1]].assemblyIndices.eqs_Fu1
                eqs_Fp1 = elements[e_assigned[n][1]].assemblyIndices.eqs_Fp1
            else # e_assigned[n][2] == 2
                eqs_Fu1 = elements[e_assigned[n][1]].assemblyIndices.eqs_Fu2
                eqs_Fp1 = elements[e_assigned[n][1]].assemblyIndices.eqs_Fp2
            end
            # Check whether this is a special node - has separate compatibility equations
            if n in specialNodesGlobalIDs       
                # Indices for new compatibility equations
                eqs_FF1 = SVector{3,Int}(i_equations+0:i_equations+2)
                eqs_FM1 = SVector{3,Int}(i_equations+3:i_equations+5)
                # Update equations' counter 
                i_equations += 6
                # Set flag
                isSpecialNode1 = true
            else
                # Set indices for compatibility equations as those already assigned
                if e_assigned[n][2] == 1
                    eqs_FF1 = elements[e_assigned[n][1]].assemblyIndices.eqs_FF1
                    eqs_FM1 = elements[e_assigned[n][1]].assemblyIndices.eqs_FM1
                else #if e_assigned[n][2] == 2
                    eqs_FF1 = elements[e_assigned[n][1]].assemblyIndices.eqs_FF2
                    eqs_FM1 = elements[e_assigned[n][1]].assemblyIndices.eqs_FM2
                end
            end
        end

        # --- Midpoint --- #
        # ----------------------------------------------------------------------

        # Indices for structural elemental states (u, p, F, M, V, Ω)
        DOF_u = SVector{3,Int}(i_states+0:i_states+2)
        DOF_p = SVector{3,Int}(i_states+3:i_states+5)
        DOF_F = SVector{3,Int}(i_states+6:i_states+8)
        DOF_M = SVector{3,Int}(i_states+9:i_states+11)
        DOF_V = SVector{3,Int}(i_states+12:i_states+14)
        DOF_Ω = SVector{3,Int}(i_states+15:i_states+17)

        # Update states' counter
        i_states += 18

        # Indices for element's generalized velocity-displacement equations
        eqs_FV = SVector{3,Int}(i_equations+0:i_equations+2)
        eqs_FΩ = SVector{3,Int}(i_equations+3:i_equations+5)

        # Update equations' counter
        i_equations += 6

        # Indices for aerodynamic states and equations
        eqs_Fχ = SVector{nAeroStates,Int}(i_equations+0:i_equations+nAeroStates-1)
        DOF_χ = SVector{nAeroStates,Int}(i_states+0:i_states+nAeroStates-1)

        # Update counters
        i_equations += nAeroStates
        i_states += nAeroStates

        # Update list with all aerodynamic DOFs
        append!(rangeAeroStates,collect(DOF_χ))

        # --- Second node --- #
        # ----------------------------------------------------------------------
            
        # Get element's second node's global ID
        n = nodesGlobalIDs[2]

        # Initialize flags for second node
        isSpecialNode2 = false
        eqsNode2Set = false

        # Check if the node's equations/states have been assigned indices
        if !assigned[n]
            # Update assignment flag and pointer to element
            assigned[n] = true
            push!(e_assigned,[e, 2])
            # Indices for equilibrium and compatibility equations of element's second node
            eqs_Fu2 = SVector{3,Int}(i_equations+0:i_equations+2)
            eqs_Fp2 = SVector{3,Int}(i_equations+3:i_equations+5)
            eqs_FF2 = SVector{3,Int}(i_equations+6:i_equations+8)
            eqs_FM2 = SVector{3,Int}(i_equations+9:i_equations+11)
            # Update equations' counter
            i_equations += 12
            # Check whether this is a special node - has nodal states
            if n in specialNodesGlobalIDs
                # Get corresponding special node
                s = findfirst(x -> n in x, specialNodesGlobalIDs)
                specialNode = specialNodes[s]
                # Indices for nodal states 
                DOF_uF = SVector{3,Int}(i_states+0:i_states+2)
                DOF_pM = SVector{3,Int}(i_states+3:i_states+5)
                # Update states' counter
                i_states += 6
                # Set flag
                isSpecialNode2 = true
                # Update special node's assembly data
                setfield!(specialNode, :assemblyData, SpecialNodeAssemblyData(specialNode, DOF_uF, DOF_pM))
            end
        else   
            # Update flag for equations having already been set
            eqsNode2Set = true
            # Set indices for equilibrium equations as those already assigned
            if e_assigned[n][2] == 1
                eqs_Fu2 = elements[e_assigned[n][1]].assemblyIndices.eqs_Fu1
                eqs_Fp2 = elements[e_assigned[n][1]].assemblyIndices.eqs_Fp1
            else #if e_assigned[n][2] == 2
                eqs_Fu2 = elements[e_assigned[n][1]].assemblyIndices.eqs_Fu2
                eqs_Fp2 = elements[e_assigned[n][1]].assemblyIndices.eqs_Fp2
            end
            # Check whether this is a special node
            if n in specialNodesGlobalIDs   
                # Indices for new compatibility equations
                eqs_FF2 = SVector{3,Int}(i_equations+0:i_equations+2)
                eqs_FM2 = SVector{3,Int}(i_equations+3:i_equations+5)
                # Update equations' counter 
                i_equations += 6
                # Set flag
                isSpecialNode2 = true
            else
                # Set indices for compatibility equations as those already assigned
                if e_assigned[n][2] == 1
                    eqs_FF2 = elements[e_assigned[n][1]].assemblyIndices.eqs_FF1
                    eqs_FM2 = elements[e_assigned[n][1]].assemblyIndices.eqs_FM1
                else #if e_assigned[n][2] == 2
                    eqs_FF2 = elements[e_assigned[n][1]].assemblyIndices.eqs_FF2
                    eqs_FM2 = elements[e_assigned[n][1]].assemblyIndices.eqs_FM2
                end
            end
        end

        # --- Update --- #
        # ----------------------------------------------------------------------
        # Update element assembly indices
        setfield!(element, :assemblyIndices, ElementAssemblyIndices{nAeroStates}(eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,eqs_Fχ,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,DOF_χ,isSpecialNode1,isSpecialNode2,eqsNode1Set,eqsNode2Set))

    end

    ## Set indices of the equations/DOFs of the hinge axis constraints
    # --------------------------------------------------------------------------
    # Initialize counter for the number of Lagrange multiplier equations in the assembly (considering all constraints)
    nλ = 0
    # Loop over hinge axis constraints
    for c in hingeAxisConstraints
        # Master and slave elements
        masterElement = c.beam.elements[c.masterLocalID]
        slaveElement = c.beam.elements[c.slaveLocalID]
        # Global ID of master and slave elements
        masterGlobalID = masterElement.assemblyData.globalID
        slaveGlobalID = slaveElement.assemblyData.globalID
        # Global indices of the moment equilibrium equations of the hinge node, associated to the master and slave elements. The element on the left keeps its equilibrium equations in the node's shared equilibrium equations (eqs_Fp2); the element on the right, whose hinged DOFs are decoupled, in the node's shared compatibility equations (eqs_FM1)
        eqs_Fp = c.masterOnLeft ?
                 vcat(masterElement.assemblyIndices.eqs_Fp2, slaveElement.assemblyIndices.eqs_FM1) :
                 vcat(masterElement.assemblyIndices.eqs_FM1, slaveElement.assemblyIndices.eqs_Fp2)
        # Global rotation DOF indices of master and slave elements
        DOF_pM = masterElement.assemblyIndices.DOF_p
        DOF_pS = slaveElement.assemblyIndices.DOF_p
        # Number of Lagrange multipliers
        N = length(c.state.λ)
        # Set range of Lagrange multiplier equations/DOFs (considering all constraints) and update counter
        λEqsRange = (nλ+1):(nλ+N)
        nλ += N
        # Set indices
        eqs_λ = i_equations - 1 .+ λEqsRange
        DOF_λ = i_states - 1 .+ λEqsRange
        # Update counters
        i_equations += length(λEqsRange)
        i_states += length(λEqsRange)
        # Update constraint assembly data
        setfield!(c, :assemblyData, HingeAxisConstraintAssemblyData(c, masterGlobalID, slaveGlobalID, eqs_Fp, DOF_pM, DOF_pS, eqs_λ, DOF_λ))
    end

    # System's order (number of equations) and nominal range
    systemOrder = i_equations - 1
    rangeNominal = 1 : systemOrder

    ## Set indices for special nodes' equations
    # --------------------------------------------------------------------------
    # Loop over special nodes
    for specialNode in specialNodes
        @unpack connectedElementsGlobalIDs,ζonElements = specialNode.assemblyData
        # Connected elements
        connectedElements = elements[connectedElementsGlobalIDs]
        # Number of connected elements
        N = length(connectedElements)
        # First connected element determines the node's equilibrium equations' indices
        eqs_Fu = ζonElements[1] == -1 ? 
                 connectedElements[1].assemblyIndices.eqs_Fu1 : 
                 connectedElements[1].assemblyIndices.eqs_Fu2
        eqs_Fp = ζonElements[1] == -1 ?
                 connectedElements[1].assemblyIndices.eqs_Fp1 :
                 connectedElements[1].assemblyIndices.eqs_Fp2                 
        # Initialize arrays for special node's compatibility equations
        eqs_FF = Vector{SVector{3,Int}}(undef, N)
        eqs_FM = Vector{SVector{3,Int}}(undef, N)
        # Loop over connected elements
        for (e, element) in enumerate(connectedElements)
            # Unpack element's assembly indices
            @unpack eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2 = element.assemblyIndices
            # Check in which side of the element the node is and set indices for nodal equations accordingly
            if ζonElements[e] == -1
                eqs_FF[e] = eqs_FF1
                eqs_FM[e] = eqs_FM1
            else # ζonElements[e] == 1
                eqs_FF[e] = eqs_FF2
                eqs_FM[e] = eqs_FM2
            end
        end
        # Update special node's assembly data
        setfield!(specialNode, :assemblyData, SpecialNodeAssemblyData(specialNode,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM))
    end
    
    return systemOrder, rangeNominal, rangeAeroStates

end


# Updates the double attachment springs data for all special nodes
function update_double_attachment_springs_data!(specialNodes)

    # Loop over special nodes
    for (n,specialNode) in enumerate(specialNodes)
        # Unpack node data
        @unpack springs = specialNode
        @unpack globalID = specialNode.assemblyData
        # Get double attachment springs on the node
        doubleAttachmentSprings = filter(spring -> spring.assemblyData.hasDoubleAttachment, springs)
        # Initialize array for double attachment spring data
        doubleAttachmentSpringsData = Vector{DoubleAttachmentSpringData}(undef,length(doubleAttachmentSprings))
        # Loop over node's springs
        for (s,spring) in enumerate(doubleAttachmentSprings)
            # Unpack spring data
            @unpack nodesGlobalIDs,nodesSpecialIDs = spring.assemblyData
            # Other node's special ID and DOF assembly indices
            otherNodeSpecialID = globalID == nodesGlobalIDs[1] ? 
                                 nodesSpecialIDs[2] : 
                                 nodesSpecialIDs[1]
            # Node of other attachment
            otherNode = specialNodes[otherNodeSpecialID] 
            # Other node's DOF assembly indices
            otherNodeDOF_uF = otherNode.assemblyData.DOF_uF
            otherNodeDOF_pM = otherNode.assemblyData.DOF_pM
            # Set spring's double attachment spring data
            doubleAttachmentSpringsData[s] = DoubleAttachmentSpringData(otherNodeSpecialID,otherNodeDOF_uF,otherNodeDOF_pM)
        end
        # Set special node's assembly data with double attachment springs data
        setfield!(specialNode, :assemblyData, SpecialNodeAssemblyData(specialNode, doubleAttachmentSpringsData))
    end

end


# Updates the assembly data of the special nodes, including trim loads and flap deflections
function update_trim_assembly_data!(beams,specialNodes,trimLoadsLinks,flapLinks,systemOrder)

    # Initialize number of independent trim loads and trim flap deflections
    nTrimLoads = 0
    nTrimδ = 0

    # Get global IDs of all nodes linked by trim loads
    allLinkedNodesGlobalIDs = reduce(vcat, [vcat(link.assemblyData.masterNodeGlobalID, link.assemblyData.slaveNodesGlobalIDs...) for link in trimLoadsLinks]; init=Int[])

    # Set independent nodal trim loads DOFs -- loop special nodes
    for specialNode in specialNodes
        # Skip nodes in trim loads links
        specialNode.assemblyData.globalID in allLinkedNodesGlobalIDs && continue
        # Initialize trim loads DoF array for the node
        DOF_T = zeros(Int, 6)
        # Loop BCs with trim loads, update nodal trim loads count and set DOF indices
        for BC in specialNode.BCs, i in eachindex(BC.deducedData.isTrim)
            if BC.deducedData.isTrim[i]
                nTrimLoads += 1
                DOF_T[i] = systemOrder + nTrimLoads
            end
        end
        # Update special node's assembly data
        setfield!(specialNode, :assemblyData, SpecialNodeAssemblyData(specialNode, SVector{6,Int}(DOF_T)))
    end

    # Set nodal trim loads DOFs for masters and respective slave nodes in trim loads links
    for link in trimLoadsLinks
        @unpack masterBC = link.masterAndSlave
        @unpack masterNodeGlobalID,slaveNodesGlobalIDs = link.assemblyData
        # Get master and slave nodes
        masterSpecialNodeID = findfirst(x -> x.assemblyData.globalID == masterNodeGlobalID, specialNodes)
        slavesSpecialNodeIDs = [findfirst(x -> x.assemblyData.globalID == slaveNodeGlobalID, specialNodes) for slaveNodeGlobalID in slaveNodesGlobalIDs]
        masterNode = specialNodes[masterSpecialNodeID]
        slaveNodes = specialNodes[slavesSpecialNodeIDs]
        # Initialize trim loads DoF array for the master and slave nodes
        DOF_T = zeros(Int, 6)
        # Update nodal trim loads count and set DOF indices through the master BC
        for i in eachindex(masterBC.deducedData.isTrim)
            if masterBC.deducedData.isTrim[i]
                nTrimLoads += 1
                DOF_T[i] = systemOrder + nTrimLoads
            end
        end
        # Update master and slave nodes' assembly data
        [setfield!(node, :assemblyData, SpecialNodeAssemblyData(node, SVector{6,Int}(DOF_T))) for node in vcat(masterNode, slaveNodes...)]
    end

    # Get IDs of all beams linked by flap deflections
    allLinkedBeamsIDs = reduce(vcat, [vcat(link.masterBeam.assemblyData.ID, [beam.assemblyData.ID for beam in link.slaveBeams]...) for link in flapLinks]; init=Int[])

    # Set trim flap deflections DOFs for independent surfaces -- loop beams
    for beam in beams
        # Unpack beam data
        @unpack ID = beam.assemblyData
        @unpack surface = beam.attachments
        # Skip if beam does not have a trim flap deflection or flap deflection is not independent
        (isnothing(surface) || !surface.assemblyData.δIsTrimVariable || ID in allLinkedBeamsIDs) && continue
        # Increment trim flap deflections count
        nTrimδ += 1
        # Update trim flap deflections DOF index
        DOF_δ = [systemOrder + nTrimLoads + nTrimδ]
        # Set DOF index for beam surface
        setfield!(surface, :assemblyData, AeroSurfaceAssemblyData(surface, DOF_δ))
    end

    # Set trim flap deflections DOFs for linked surfaces -- loop flap links
    for link in flapLinks
        # Unpack beam data
        masterSurface = link.masterBeam.attachments.surface
         # Skip if master beam does not have a trim flap deflection
        (isnothing(masterSurface) || !masterSurface.assemblyData.δIsTrimVariable) && continue
        # Increment trim flap deflections count
        nTrimδ += 1
        # Update trim flap deflections DOF index
        DOF_δ = [systemOrder + nTrimLoads + nTrimδ]
        # Set DOF index for master surface
        setfield!(masterSurface, :assemblyData, AeroSurfaceAssemblyData(masterSurface, DOF_δ))
        # Set DOF index for slave surfaces
        [setfield!(slaveSurface, :assemblyData, AeroSurfaceAssemblyData(slaveSurface, DOF_δ)) for slaveSurface in (beam.attachments.surface for beam in link.slaveBeams)]
    end
    
    # Total number of trim variables
    nTrimVariables = nTrimLoads + nTrimδ

    # Set ranges
    rangeTrimLoads = systemOrder + 1 : systemOrder + nTrimLoads
    rangeTrimδ = systemOrder + nTrimLoads + 1 : systemOrder + nTrimLoads + nTrimδ

    return nTrimVariables, nTrimLoads, nTrimδ, rangeTrimLoads, rangeTrimδ

end


# Updates the initial condition states on all elements and nodes of the assembly
function update_initial_states!(beams,elements,specialNodes,kinematicsBasisA,atmosphere)

    # Loop over elements
    for element in elements
        
        # Unpack
        @unpack Δℓ,x1,x1_n1,x1_n2,k,R0,R0T,R0T_n1,R0T_n2 = element.geometry
        @unpack parentBeamID,r = element.assemblyData
        @unpack C,I = element.sectional
        @unpack u,p,F,M,V,Ω,γ,κ,P,H,χ = element.state.midpoint
        @unpack udot,pdot,Vdot,Ωdot = element.state.rates
        @unpack u1,u2,p1,p2,u1_b,u2_b,p1_b,p2_b,F1,F2,M1,M2,θ = element.state.nodal
        @unpack R,HT,RR0,RR0T = element.stateCache.rotationMatrices
        @unpack v,ω,vdot,ωdot = element.stateCache.kinematicsBasisb
        @unpack v_A,ω_A,vdot_A,ωdot_A,R_A,R_AT = kinematicsBasisA.state
        @unpack u0_of_x1,p0_of_x1,udot0_of_x1,pdot0_of_x1,uprime0_of_x1,pprime0_of_x1 = beams[parentBeamID].initialConditions

        # Generalized velocities of basis b at element's midpoint, resolved in basis A
        v .= R_AT * (v_A + cross(ω_A,R_A*r))
        ω .= R_AT * ω_A

        # Generalized accelerations of basis b at element's midpoint, resolved in basis A
        vdot .= R_AT * (vdot_A + cross(ωdot_A,R_A*r) + cross(ω_A,cross(ω_A,R_A*r)) - cross(ω_A,R_A*v_A))
        ωdot .= R_AT * ωdot_A

        # Generalized displacements, resolved in basis A
        u .= R0 * u0_of_x1(x1)
        p .= R0 * p0_of_x1(x1)

        # Derivatives of generalized displacements, resolved in basis A
        u_prime = R0 * uprime0_of_x1(x1)
        p_prime = R0 * pprime0_of_x1(x1)

        # Rotation and tangent tensors, resolved in basis A
        R .= rotation_tensor_WM(p)
        HT .= HT_WM(p)
        RR0 .= R*R0
        RR0T .= RR0'

        # Sectional strains, resolved in basis b
        γ .= (R-I3) * a1 + RR0T * u_prime
        κ .= R0T * HT * p_prime
        sectionalStrains = [γ; κ]

        # Sectional forces, resolved in basis B
        sectionalForces = inv(C) * sectionalStrains
        F .= sectionalForces[1:3]
        M .= sectionalForces[4:6]

        # Generalized displacements' rates, resolved in basis A
        udot .= R0 * udot0_of_x1(x1)
        pdot .= R0 * pdot0_of_x1(x1)

        # Generalized sectional velocities, resolved in basis B
        V .= RR0T * (v + cross(ω,u) + udot)
        Ω .= R0T * (HT * pdot + R' * ω)
        sectionalVelocities = [V; Ω]
        
        # Sectional linear and angular momenta, resolved in basis B 
        sectionalMomenta = I * sectionalVelocities
        P .= sectionalMomenta[1:3]
        H .= sectionalMomenta[4:6]

        # Midpoint deformed curvature vector, resolved in basis B
        K = k + κ

        # Midpoint generalized forces' derivatives with respect to x1 (neglecting the unknown momenta rates)
        F_prime = cross(Ω,P) - cross(K,F)
        M_prime = cross(Ω,H) + cross(V,P) - cross(a1+γ,F) - cross(K,M)

        # Nodal displacements and rotations, resolved in basis A
        u1 .= u0_of_x1(x1_n1)
        u2 .= u0_of_x1(x1_n2)
        p1 .= p0_of_x1(x1_n1)
        p2 .= p0_of_x1(x1_n2)

        # Nodal displacements and rotation parameters' vectors, resolved in basis b
        u1_b .= R0T_n1 * u1
        u2_b .= R0T_n2 * u2
        p1_b .= R0T_n1 * p1
        p2_b .= R0T_n2 * p2

        # Nodal rotation angles
        θ1 = rotation_angle(p1)
        θ2 = rotation_angle(p2)
        θ .= [θ1; θ2]

        # Nodal forces and moments (neglecting distributed loads' resultants)
        F1 .= F - (Δℓ/2 * F_prime)
        F2 .= F + (Δℓ/2 * F_prime)
        M1 .= M - (Δℓ/2 * M_prime)
        M2 .= M + (Δℓ/2 * M_prime)

        # Aerodynamic states
        χ .= initial_aero_states!(atmosphere,element,χ,V,Ω,Vdot,Ωdot)

    end

    # Loop over special nodes
    for specialNode in specialNodes
        @unpack connectedElementsGlobalIDs,ζonElements = specialNode.assemblyData
        @unpack u,p,F,M = specialNode.state
        connectedElements = elements[connectedElementsGlobalIDs]
        # Update special node's state with the first connected element's nodal states
        u .= ζonElements[1] == -1 ?
             connectedElements[1].state.nodal.u1 :
             connectedElements[1].state.nodal.u2
        p .= ζonElements[1] == -1 ?
             connectedElements[1].state.nodal.p1 :
             connectedElements[1].state.nodal.p2
        F .= ζonElements[1] == -1 ?
             connectedElements[1].state.nodal.F1 :
             connectedElements[1].state.nodal.F2
        M .= ζonElements[1] == -1 ?
             connectedElements[1].state.nodal.M1 :
             connectedElements[1].state.nodal.M2
    end

    # Set initial state of the model
    initialState = InitialModelState(
        [copy(element.state) for element in elements],
        [copy(specialNode.state) for specialNode in specialNodes],
    )

    return initialState

end


# Helper function to update the model's assembly data
function update_model_assembly!(model)

    @unpack beams,connections,firstNodePosition,trimLoadsLinks,flapLinks,hingeAxisConstraints,BCs,kinematicsBasisA,atmosphere,gust = model

    # Get new assembly data
    assemblyData,BCs = get_model_assembly_data!(beams,connections,firstNodePosition,trimLoadsLinks,flapLinks,hingeAxisConstraints,BCs,gust)

    # Update initial states
    initialState = update_initial_states!(beams,assemblyData.elements,assemblyData.specialNodes,kinematicsBasisA,atmosphere)

    # Set on model
    setfield!(model, :assemblyData, assemblyData)
    setfield!(model, :BCs, BCs)
    setfield!(model, :initialState, initialState)

end


# Resets the model's elements' and special nodes' states to their initial values
function reset_model_state!(model)

    # Unpack
    @unpack elements,specialNodes = model.assemblyData

    # Loop over elements and reset their states
    for (element,initialState) in zip(elements,model.initialState.elementsStates)
        setfield!(element, :state, copy(initialState))
    end

    # Loop over special nodes and reset their states
    for (specialNode,initialState) in zip(specialNodes,model.initialState.specialNodesStates)
        setfield!(specialNode, :state, copy(initialState))
    end

end