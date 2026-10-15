## BEAM COMPOSITE TYPE
# ------------------------------------------------------------------------------
const _BEAM_ID_COUNTER = Ref(0)
next_beam_construction_id!() = (_BEAM_ID_COUNTER[] += 1)

# Beam geometry composite type
# Fields
# - `ℓ`: (arc)length of the beam
# - `r`: position vector of the beam's endpoint relative to its starting point, resolved in the beam's undeformed basis b(0)
# - `firstNodeOffset`: position offset of the beam's first node relative to the its connection node (if any)
# - `k`: undeformed beam's curvature vector per unit length
# - `p0`: Wiener-Milenkovic rotation parameters from basis A to basis b(0)
# - `R0`: rotation tensor associated with p0
#
struct BeamGeometry

    ℓ::Float64
    r::SVector{3,Float64}
    firstNodeOffset::SVector{3,Float64}
    k::SVector{3,Float64}
    p0::SVector{3,Float64}
    R0::SMatrix{3,3,Float64,9}

end
    

# Beam sectional data composite type
# Fields
# - `S`: array of sectional stiffness matrices
# - `I`: array of sectional inertia matrices
#
struct BeamSectionalData

    S::Vector{SMatrix{6,6,Float64,36}}
    I::Vector{SMatrix{6,6,Float64,36}}

end


# Beam discretization composite type
# Fields
# - `nElements`: number of elements
# - `normalizedNodalPositions`: normalized nodal positions of beam elements
#
struct BeamDiscretization

    nElements::Int
    normalizedNodalPositions::Vector{Float64}

end


# Beam hinges' data composite type
# Fields
# - `hingedNodes`: nodes with a hinge
# - `hingedNodesDoF`: respective hinged degrees-of-freedom
#
struct BeamHingesData

    hingedNodes::Vector{Int}
    hingedNodesDoF::Vector{SVector{3,Bool}}

end


# Beam attachments composite type
# Fields
# - `springs`: array of attached springs
# - `pointInertias`: array of attached point inertias
# - `surface`: attached aerodynamic surface
#
mutable struct BeamAttachments

    springs::Vector{Spring}
    pointInertias::Vector{PointInertia}
    surface::Union{Nothing,AeroSurface}
    
end
BeamAttachments(surface) = BeamAttachments(
    Spring[], PointInertia[], surface
    )


# Beam initial conditions composite type
# Fields
# - `hasInitialConditions`: flag for having initial conditions
# - `u0_of_x1`: initial displacement (resolved in the undeformed beam basis, b) of the beam as a function of its arclength coordinate (x1)
# - `p0_of_x1`: initial rotation parameters (resolved in the undeformed beam basis, b) of the beam as a function of its arclength coordinate (x1)
# - `udot0_of_x1`: initial displacement's rates (resolved in the undeformed beam basis, b) of the beam as a function of its arclength coordinate (x1)
# - `pdot0_of_x1`: initial rotation parameters' rates (resolved in the undeformed beam basis, b) of the beam as a function of its arclength coordinate (x1)
# - `uprime0_of_x1`: derivative of initial displacement with respect to arclength coordinate (x1)
# - `pprime0_of_x1`: derivative of initial rotation parameters with respect to arclength coordinate (x1)
#
mutable struct BeamInitialConditions

    hasInitialConditions::Bool

    u0_of_x1::Function
    p0_of_x1::Function

    udot0_of_x1::Function
    pdot0_of_x1::Function

    uprime0_of_x1::Function
    pprime0_of_x1::Function

end
BeamInitialConditions() = BeamInitialConditions(
    false,
    ntuple(_ -> x1 -> z3, 2)...,
    ntuple(_ -> x1 -> z3, 2)...,
    ntuple(_ -> x1 -> z3, 2)...,
    )


# Beam distributed loads composite type
# Fields
# - `f_A_of_x1t`: distributed forces in basis A as a function of arclength coordinate (x1) and time (t)
# - `m_A_of_x1t`: distributed moments in basis A as a function of arclength coordinate (x1) and time (t)
# - `f_b_of_x1t`: distributed forces in basis b as a function of arclength coordinate (x1) and time (t)
# - `m_b_of_x1t`: distributed moments in basis b as a function of arclength coordinate (x1) and time (t)
# - `ff_A_of_x1t`: follower forces in basis A as a function of arclength coordinate (x1) and time (t)
# - `mf_A_of_x1t`: follower moments in basis A as a function of arclength coordinate (x1) and time (t)
# - `ff_b_of_x1t`: follower forces in basis b as a function of arclength coordinate (x1) and time (t)
# - `mf_b_of_x1t`: follower moments in basis b as a function of arclength coordinate (x1) and time (t)
#
mutable struct BeamDistributedLoads

    f_A_of_x1t::Union{Nothing,<:Function}
    m_A_of_x1t::Union{Nothing,<:Function}
    f_b_of_x1t::Union{Nothing,<:Function}
    m_b_of_x1t::Union{Nothing,<:Function}
    ff_A_of_x1t::Union{Nothing,<:Function}
    mf_A_of_x1t::Union{Nothing,<:Function}
    ff_b_of_x1t::Union{Nothing,<:Function}
    mf_b_of_x1t::Union{Nothing,<:Function}

end
BeamDistributedLoads() = BeamDistributedLoads(
    ntuple(_ -> nothing, 8)...,
    )
    

# Beam assembly data composite type
# Fields
# - `ID`: global ID of the beam in the assembly
# - `elementsGlobalIDs`: range of element global IDs belonging to the beam
# - `nodesGlobalIDs`: range of node global IDs belonging to the beam
# - `r_n`: position vectors of the nodes relative to the assembly's origin (basis A)
# - `connectedBeams`: array of beams to which this beam is connected (a non-recursive property)
# - `connectedNodesThis`: nodes of this beam which are connected to other beams' nodes
# - `connectedNodesOther`: respective nodes of the other beams
#
struct BeamAssemblyData

    ID::Int

    elementsGlobalIDs::Vector{Int}
    nodesGlobalIDs::Vector{Int}

    r_n::Vector{SVector{3,Float64}}

    connectedBeamsIDs::Vector{Int}
    connectedNodesThis::Vector{Int}
    connectedNodesOther::Vector{Int}

end
BeamAssemblyData() = BeamAssemblyData(
    0, 
    Vector{Int}(), Vector{Int}(), 
    Vector{SVector{3,Float64}}(),
    Vector{Int}(), Vector{Int}(), Vector{Int}()
    )


# Beam composite type
# Fields
# - `name`: name of the beam
# - `constructionID`: construction ID of the beam
# - `geometry`: beam geometry
# - `sectional`: beam sectional properties
# - `discretization`: beam discretization
# - `initialConditions`: beam initial conditions
# - `distLoads`: beam distributed loads
# - `attachments`: beam attachments (springs, point inertias, aerodynamic surface)
# - `elements`: array of beam elements
# - `assemblyData`: assembly data of the beam
#
@with_kw mutable struct Beam

    # Name and construction ID
    name::String
    constructionID::Int

    # Geometry
    geometry::BeamGeometry

    # Sectional properties
    sectional::BeamSectionalData

    # Discretization
    discretization::BeamDiscretization

    # Hinges data
    hinges::BeamHingesData

    # Attachments
    attachments::BeamAttachments

    # Elements
    elements::Vector{Element}

    # Initial conditions
    initialConditions::BeamInitialConditions = BeamInitialConditions()

    # Distributed loads
    distLoads::BeamDistributedLoads = BeamDistributedLoads()

    # Assembly data (to be set by Model)
    assemblyData::BeamAssemblyData = BeamAssemblyData()

end
export Beam


## PUBLIC API
# ------------------------------------------------------------------------------
"""
    create_Beam(; kwargs...)

Beam constructor

# Keyword Arguments
- `name::String`: name of the beam
- `length::Union{Missing,<:Real}`: (arc)length of the beam
- `r::Union{Missing,AbstractVector{<:Real}}`: position vector of the beam's endpoint relative to the assembly's origin (basis A)
- `k::AbstractVector{<:Real}`: undeformed beam's curvature vector per unit length
- `p0::AbstractVector{<:Real}`: Wiener-Milenkovic rotation parameters from basis A to basis b(0)
- `S::AbstractVector{<:AbstractMatrix{<:Real}}`: array of sectional stiffness matrices
- `I::AbstractVector{<:AbstractMatrix{<:Real}}`: array of sectional inertia matrices
- `nElements::Int`: number of elements for discretization
- `normalizedNodalPositions::AbstractVector{<:Real}`: normalized nodal positions of beam elements
- `hingedNodes::AbstractVector{Int}`: nodes with a hinge
- `hingedNodesDoF::AbstractVector{<:AbstractVector{Bool}}`: respective hinged degrees-of-freedom
- `surface::Union{Nothing,AeroSurface}`: attached aerodynamic surface
"""
function create_Beam(; 
    name::String="",
    length::Union{Missing,<:Real}=missing,
    r::Union{Missing,AbstractVector{<:Real}}=missing,
    firstNodeOffset::AbstractVector{<:Real}=zeros(3),
    k::AbstractVector{<:Real}=zeros(3),
    p0::AbstractVector{<:Real}=zeros(3),
    S::AbstractVector{<:AbstractMatrix{<:Real}},
    I::Union{Missing,AbstractVector{<:AbstractMatrix{<:Real}}}=missing,
    nElements::Int,
    normalizedNodalPositions::AbstractVector{<:Real}=Vector{Float64}(),
    hingedNodes::AbstractVector{Int}=Vector{Int}(),
    hingedNodesDoF::AbstractVector{<:AbstractVector{Bool}}=Vector{Vector{Bool}}(),
    surface::Union{Nothing,AeroSurface}=nothing
    )

    # Validate inputs and get beam fields
    geometry,discretization,sectional,hinges,attachments = validate_and_set_beam(length,r,firstNodeOffset,k,p0,S,I,nElements,normalizedNodalPositions,hingedNodes,hingedNodesDoF,surface)

    # Create beam elements
    elements = create_beam_elements(geometry,sectional,discretization,hinges,attachments)

    return Beam(
        name=name,
        constructionID=next_beam_construction_id!(),
        geometry=geometry,
        sectional=sectional,
        discretization=discretization,
        hinges=hinges,
        attachments=attachments,
        elements=elements,
    )

end
export create_Beam


"""
    update_beam_orientation!(; beam, p0)

Updates the orientation of the beam

# Keyword arguments
- `beam::Beam`: reference beam
- `p0::AbstractVector{<:Real}`: Wiener-Milenkovic rotation parameters from basis A to basis b(0)
"""
function update_beam_orientation!(; beam::Beam,p0::AbstractVector{<:Real})

    # Validate and get new geometry
    geometry = validate_and_set_geometry(beam.geometry.ℓ,beam.geometry.r,beam.geometry.firstNodeOffset,beam.geometry.k,p0)

    # Update beam geometry field
    setfield!(beam, :geometry, geometry)

    # Update elements' geometry
    for element in beam.elements
       setfield!(element, :geometry, element_geometry(element.assemblyData.nodesLocalIDs,beam.geometry.ℓ,beam.geometry.k,beam.geometry.R0,beam.discretization.normalizedNodalPositions,beam.hinges.hingedNodes,beam.hinges.hingedNodesDoF))
    end

end
export update_beam_orientation!


"""
    add_distributed_loads_to_beam!(; beam, loads)

Adds loads to the beam

# Keyword arguments
- `beam::Beam`: reference beam
- `loads::Dict{Symbol,<:Function}`: dictionary whose keys are load types and values are the respective load functions
"""
function add_distributed_loads_to_beam!(; beam::Beam, loads::Dict{Symbol,<:Function})
    
    # Loop loads
    for (load,fun) in loads

        # Check inputs
        @assert load in [:f_A_of_x1t,:m_A_of_x1t,:f_b_of_x1t,:m_b_of_x1t,:ff_A_of_x1t,:mf_A_of_x1t,:ff_b_of_x1t,:mf_b_of_x1t]
        @assert applicable(fun, 0, 0) "function for $load must be callable with two arguments (x1,t)"
        @assert isa(fun(0,0),AbstractVector{<:Real}) "$load must return a vector of real numbers when called with two arguments (x1,t)"
        @assert length(fun(0,0)) == 3 "$load must return a vector of length 3 when called with two arguments (x1,t)"

        # Update load on beam
        setfield!(beam.distLoads, load, fun)
        
        # Update load on elements
        for element in beam.elements
            update_element_distributed_loads!(element,load,fun)
        end

    end

end
export add_distributed_loads_to_beam!


"""
    add_initial_conditions_to_beam!(; beam,conditions)

Adds initial generalized displacements and velocities to the beam. Any combination of them can be specified: the components that are not specified (null) are obtained from the specified ones through the governing equations, so that the initial states of the dynamic problem are always consistent. In particular, the initial states are such that:
- The specified displacements and rotations are satisfied, and the unspecified ones are those of the static equilibrium (which makes the initial rotations compatible with the initial displacements, and vice-versa)
- The specified velocities are satisfied, and the unspecified ones follow the specified ones as in a quasi-static motion (e.g., the angular velocities that accompany specified displacements' rates, and vice-versa)
- The initial accelerations satisfy the equations of motion

# Keyword arguments
- `beam::Beam`: reference beam
- `conditions::Dict{Symbol,<:Function}`: dictionary whose keys are initial condition types (`:u0_of_x1`, `:p0_of_x1`, `:udot0_of_x1` and `:pdot0_of_x1`) and values are the respective functions over arclength. These return the initial displacements, rotation parameters and their rates, respectively, resolved in the undeformed beam basis (b)
"""
function add_initial_conditions_to_beam!(; beam::Beam, conditions::Dict{Symbol,<:Function})
    
    # Loop initial condition types and respective functions
    for (type,fun) in conditions

        # Check inputs
        @assert type in [:u0_of_x1,:p0_of_x1,:udot0_of_x1,:pdot0_of_x1]
        @assert applicable(fun, 0) "function for $type must be callable with one argument (x1)"
        @assert isa(fun(0),AbstractVector{<:Real}) "$type must return a vector of real numbers when called with a single argument (x1)"
        @assert length(fun(0)) == 3 "$type must return a vector of length 3 when called with a single argument (x1)"

        # Set initial condition on beam
        setfield!(beam.initialConditions, type, fun)

        # The update of elemental states occurs on the model creation/update (because V and Ω depend on the model's v_A and ω_A)
    end

    # Get derivatives of initial generalized displacements with respect to arclength coordinate
    uprime0_of_x1, pprime0_of_x1 = initial_displacements_derivatives(beam.initialConditions.u0_of_x1,beam.initialConditions.p0_of_x1)

    # Set on beam
    setfield!(beam.initialConditions, :uprime0_of_x1, uprime0_of_x1)
    setfield!(beam.initialConditions, :pprime0_of_x1, pprime0_of_x1)
    setfield!(beam.initialConditions, :hasInitialConditions, true)

end
export add_initial_conditions_to_beam!


"""
    add_point_inertia_to_beam!(; beam,elementID,pointInertia)

Adds a point inertia to a particular element of the beam

# Keyword arguments
- `beam::Beam`: reference beam
- `elementID::Int`: local ID (on the beam) of the element to which the point inertia is attached
- `pointInertia::PointInertia`: point inertia
"""
function add_point_inertia_to_beam!(; beam::Beam,elementID::Int,pointInertia::PointInertia)

    # Check that the beam has the element to which the point inertia was assigned
    @assert elementID <= beam.discretization.nElements "point inertia element ID out of range for beam with $(beam.discretization.nElements) elements"

    # Add point inertia to beam attachments
    push!(beam.attachments.pointInertias, pointInertia)

    # Add point inertia to element and update its sectional data
    add_point_inertia_to_element!(beam.elements[elementID],pointInertia)

end
export add_point_inertia_to_beam!


"""
    add_spring_to_beam!(; beam, springs, nodeID, p, basis)

Adds a simply-attached spring to a beam

# Keyword arguments
- `beam::Beam`: reference beam
- `spring::Spring`: spring
- `nodeID::Int`: local ID (on the beam) of the node to which the spring is attached
- `Δu0::Real`: translational pre-deformation
- `Δp0::Real`: rotational pre-deformation
- `p::AbstractVector{<:Real}`: Euler 3-2-1 rotation parameters from basis b of the attachment node to the spring's basis s
- `basis::Symbol`: basis in which the spring's attachment is defined, either :b or :A
"""
function add_spring_to_beam!(; beam::Beam,spring::Spring,nodeID::Int,Δu0::Real=0,Δp0::Real=0,p::AbstractVector{<:Real}=zeros(3),basis::Symbol=:b)

    # Validate
    @assert nodeID > 0 "nodeID must be positive"
    @assert nodeID <= beam.discretization.nElements+1 "nodeID out of range for beam with $(beam.discretization.nElements+1) nodes"
    @assert length(p) == 3 "p must be a 3-element vector"
    @assert basis == :b || basis == :A "basis must be either :b or :A"
    Δu0 = Float64(Δu0)
    Δp0 = Float64(Δp0)

    # Update spring assembly data
    update_spring_assembly_data!([beam],spring,[nodeID],Δu0,Δp0,p,basis)

    # Add spring to beam attachments
    push!(beam.attachments.springs, spring)

end
export add_spring_to_beam!


"""
    add_spring_to_beams!(; beams, spring, nodesLocalIDs, p, Δu0, Δp0, basis)

Adds a doubly-attached spring to the beams

# Keyword arguments
- `beams::Vector{Beam}`: reference beams
- `spring::Spring`: doubly-attached spring
- `nodesLocalIDs`: local IDs (on the respective beams) of the nodes to which the spring is attached
- `Δu0::Real`: translational pre-deformation
- `Δp0::Real`: rotational pre-deformation
- `p::AbstractVector{<:Real}`: Euler 3-2-1 rotation parameters from basis b of the attachment node to the spring's basis s
- `basis::Symbol`: basis in which the spring's attachment is defined, either :b or :A
"""
function add_spring_to_beams!(; beams::Vector{Beam},spring::Spring,nodesLocalIDs::AbstractVector{Int},Δu0::Real=0,Δp0::Real=0,p::AbstractVector{<:Real}=zeros(3),basis::Symbol=:b)

    # Validate
    @assert length(beams) == 2 "spring must be attached to two beams"
    @assert length(nodesLocalIDs) == 2 "spring must have two nodesLocalIDs, one for each of the two beams"
    @assert all(nodesLocalIDs .> 0) "nodesLocalIDs must be positive"
    for i in 1:2
        @assert nodesLocalIDs[i] <= beams[i].discretization.nElements+1 "nodesLocalIDs out of range for beam $(i) with $(beams[i].discretization.nElements+1) nodes"
    end
    @assert length(p) == 3 "p must be a 3-element vector"
    @assert beams[1] === beams[2] && length(unique(nodesLocalIDs)) == 2 && iszero(beams[1].geometry.k[2]) && iszero(beams[1].geometry.k[3]) "currently only supports attaching a spring to the same uncurved beam, at two different nodes"
    @assert basis == :b || basis == :A "basis must be either :b or :A"
    Δu0 = Float64(Δu0)
    Δp0 = Float64(Δp0)

    # Update spring assembly data
    update_spring_assembly_data!(beams,spring,nodesLocalIDs,Δu0,Δp0,p,basis)

    # Loop beams
    for i in 1:2
        # Add spring to beam attachments
        push!(beams[i].attachments.springs, spring)
    end

end
export add_spring_to_beams!


"""
    add_surface_to_beam!(; beam, surface)

Adds an aerodynamic surface to the beam

# Keyword arguments
- `beam::Beam`: reference beam
- `surface::AeroSurface`: aerodynamic surface
"""
function add_surface_to_beam!(; beam::Beam, surface::AeroSurface)

    # Check that there is no surface already attached to the beam
    @assert isnothing(beam.attachments.surface) "beam already has an aerodynamic surface attached, remove it first before adding a new one"

    # Update aerodynamic surface assembly data and attach to beam
    setfield!(beam.attachments, :surface, update_aero_surface!(beam.geometry.ℓ,surface))

    # Update aerodynamic surface on elements and reset the elements in the beam elements array
    for element in beam.elements
        setfield!(element, :aero, AeroProperties(surface,element.geometry))
    end

end
export add_surface_to_beam!


"""
    remove_all_loads_from_beams!(; beams)

Removes all loads from the beams

# Keyword arguments
- `beams::AbstractVector{Beam}`: array of beams
"""
function remove_all_loads_from_beams!(; beams::AbstractVector{Beam})

    # Remove from beams and their elements
    for beam in beams
        for field in [:f_A_of_x1t,:m_A_of_x1t,:f_b_of_x1t,:m_b_of_x1t,:ff_A_of_x1t,:mf_A_of_x1t,:ff_b_of_x1t,:mf_b_of_x1t]
            setfield!(beam.distLoads, field, nothing)
        end
        for element in beam.elements
            for field in [:f_A_of_ζt,:m_A_of_ζt,:f_b_of_ζt,:m_b_of_ζt,:ff_A_of_ζt,:mf_A_of_ζt,:ff_b_of_ζt,:mf_b_of_ζt]
                setfield!(element.distLoads, field, (ζ,t) -> z3)
            end
            for field in [:has_f_A,:has_m_A,:has_f_b,:has_m_b,:has_ff_A,:has_mf_A,:has_ff_b,:has_mf_b]
                setfield!(element.distLoads, field, false)
            end
        end
    end

end
export remove_all_loads_from_beams!


"""
    remove_all_initial_conditions_from_beams!(; beams)

Removes all initial conditions from the beams

# Keyword arguments
- `beams::AbstractVector{Beam}`: array of beams
"""
function remove_all_initial_conditions_from_beams!(; beams::AbstractVector{Beam})

    for beam in beams
        for field in [:u0_of_x1,:p0_of_x1,:udot0_of_x1,:pdot0_of_x1,:uprime0_of_x1,:pprime0_of_x1]
            setfield!(beam.initialConditions, field, x1 -> z3)
        end
        setfield!(beam.initialConditions, :hasInitialConditions, false)
    end

end


"""
    remove_all_point_inertias_from_beams!(; beams)

Removes all point inertias attached to the beams

# Keyword arguments
- `beams::AbstractVector{Beam}`: array of beams
"""
function remove_all_point_inertias_from_beams!(; beams::AbstractVector{Beam})

    for beam in beams
        setfield!(beam.attachments, :pointInertias, PointInertia[])
    end

end
export remove_all_point_inertias_from_beams!


"""
    remove_all_springs_from_beams!(; beams)

Removes all springs attached to the beams

# Keyword arguments
- `beams::AbstractVector{Beam}`: array of beams
"""
function remove_all_springs_from_beams!(; beams::AbstractVector{Beam})

    for beam in beams
        setfield!(beam.attachments, :springs, Spring[])
    end

end
export remove_all_springs_from_beams!


"""
    remove_surface_from_beams!(; beams)

Removes the aerodynamic surface attached to the beams

# Keyword arguments
- `beams::AbstractVector{Beam}`: array of beams
"""
function remove_surface_from_beams!(; beams::AbstractVector{Beam})

    for beam in beams
        # Check that there is a surface attached to the beam
        @assert !isnothing(beam.attachments.surface) "beam has no aerodynamic surface attached"

        # Remove aerodynamic surface from beam
        setfield!(beam.attachments, :surface, nothing)

        # Remove aerodynamic surface from elements
        for element in beam.elements
            setfield!(element, :aero, nothing)
        end
    end

end
export remove_surface_from_beams!


## PRIVATE API
# -----------------------------------------------------------------------------
# Validates the beam inputs and outputs beam fields
function validate_and_set_beam(ℓ,r,firstNodeOffset,k,p0,S,I,nElements,normalizedNodalPositions,hingedNodes,hingedNodesDoF,surface)

    # Beam creation pre-validation checks
    @assert ismissing(ℓ) != ismissing(r) "only one of either beam length or r must be specified"

    # Validate and set geometry
    geometry = validate_and_set_geometry(ℓ,r,firstNodeOffset,k,p0)

    # Validate and set discretization
    discretization = validate_and_set_discretization(nElements,normalizedNodalPositions)

    # Validate and set sectional data
    sectional = validate_and_set_sectional_data(S,I,nElements)

    # Validate and set hinges data
    hinges = validate_and_set_hinges(hingedNodes,hingedNodesDoF,nElements)

    # Update beam-dependent properties of attached aerodynamic surface
    update_aero_surface!(ℓ,surface)

    # Initialize beam attachments
    attachments = BeamAttachments(surface)

    return geometry,discretization,sectional,hinges,attachments

end


# Checks that the beam length and curvature vector are valid
function validate_and_set_geometry(ℓ,r,firstNodeOffset,k,p0)

    # Validate
    @assert ismissing(ℓ) || ℓ > 0 "beam length must be greater than zero"
    @assert ismissing(r) || length(r) == 3 "r must be a 3-element vector"
    @assert length(firstNodeOffset) == 3 "firstNodeOffset must be a 3-element vector"
    @assert length(k) == 3 "curvature vector k must be a 3-element vector"
    @assert length(p0) == 3 "p0 must have length 3"

    # Get rotation tensor from basis A to basis b(0)
    R0 = rotation_tensor_WM(p0)

    # Compute beam length from endpoint vector, or endpoint vector from beam length
    if ismissing(ℓ)
        ℓ = beam_length(r,k)
    else
        r = position_vector(R0,k,ℓ)
    end    

    # Set correct types
    ℓ = Float64(ℓ)
    r = SVector{3,Float64}(r)
    firstNodeOffset = SVector{3,Float64}(firstNodeOffset)
    k = SVector{3,Float64}(k)
    p0 = SVector{3,Float64}(p0)
    R0 = SMatrix{3,3,Float64,9}(R0)

    # Set geometry
    geometry = BeamGeometry(ℓ,r,firstNodeOffset,k,p0,R0)

    return geometry

end


# Checks that the discretization (number of elements and normalized nodal positions) is valid
function validate_and_set_discretization(nElements,normalizedNodalPositions)

    # Validate
    @assert nElements > 0 "number of elements must be greater than zero"
    if isempty(normalizedNodalPositions)
        normalizedNodalPositions = LinRange(0.0,1.0,nElements+1)
    else
        @assert iszero(first(normalizedNodalPositions)) "first normalized nodal position must be 0"
        @assert last(normalizedNodalPositions) == 1.0 "last normalized nodal position must be 1.0"
        @assert length(normalizedNodalPositions) == nElements+1 "number of normalized nodal positions must be equal to number of elements + 1"
        @assert issorted(normalizedNodalPositions) "normalized nodal positions must be sorted in ascending order"
    end

    # Set correct types
    normalizedNodalPositions = Float64.(normalizedNodalPositions)

    # Set discretization
    discretization = BeamDiscretization(nElements,normalizedNodalPositions)

    return discretization

end


# Checks that sectional matrices are input as a single one for the whole beam or are input in a per element basis 
function validate_and_set_sectional_data(S,I,nElements)

    # Sectional stiffness matrix
    @assert (length(S)==1 || length(S)==nElements) "input either one stiffness matrix for the entire beam, or one for each element"
    for Si in S
        @assert size(Si) == (6,6) "stiffness matrices must be of size (6,6)"
        @assert all([Si[j,j] > 0 for j in 1:6]) "diagonal elements of stiffness matrices must be positive"
    end

    # Set a default sectional inertia matrix if not provided (not needed for steady problems)
    if ismissing(I)
        I = [inertia_matrix(ρA=1,ρIy=1,ρIz=1)]
    end

    # Sectional inertia matrix
    @assert (length(I)==1 || length(I)==nElements) "input either one inertia matrix for the entire beam, or one for each element"
    for Ii in I
        @assert size(Ii) == (6,6) "inertia matrices must be of size (6,6)"
        @assert all([Ii[j,j] >= 0 for j in 1:6]) "diagonal elements of inertia matrices must be greater than or equal to zero"
    end

    # Set correct types
    S = SMatrix{6,6,Float64,36}.(S)
    I = SMatrix{6,6,Float64,36}.(I)

    # Set sectional data
    sectional = BeamSectionalData(S,I)

    return sectional

end


# Validates the hinged nodes and DoFs
function validate_and_set_hinges(hingedNodes,hingedNodesDoF,nElements)

    # Validate and set correct types
    if !isempty(hingedNodes)
        @assert length(hingedNodes) == length(hingedNodesDoF) "hingedNodes and hingedNodesDoF arrays must have the same length"
        @assert all(x -> 0 < x <= nElements+1, hingedNodes) "all hinged nodes must be greater than zero and less than or equal to number of elements + 1"
        @assert !any(hingedNodes[i]+1 == hingedNodes[i+1] for i in 1:length(hingedNodes)-1) "cannot have two consecutive nodes hinged"
        for (i,hingedNodeDoF) in enumerate(hingedNodesDoF)
            @assert length(hingedNodeDoF) == 3 "each vector of hingedNodesDoF must have a length of 3"
        end
        hingedNodes = collect(Int, hingedNodes)
        hingedNodesDoF = SVector{3,Bool}.(hingedNodesDoF)
    end

    # Set hinges data
    hinges = BeamHingesData(hingedNodes,hingedNodesDoF)

    return hinges

end


# Updates the beam-dependent properties of the attached aerodynamic surface 
function update_aero_surface!(ℓ,surface)

    # Skip if there is no surface
    isnothing(surface) && return surface

    # Unpack surface data
    @unpack c,hasSymmetricCounterpart = surface.geometry

    # Surface planform area
    area = c isa Real ?
           c * ℓ :
           gauss_legendre(c, 0, ℓ)

    # Aspect ratio
    AR = hasSymmetricCounterpart ? 
         2*ℓ^2/area :
         ℓ^2/area

    # Update aerodynamic surface assembly data
    setfield!(surface, :assemblyData, AeroSurfaceAssemblyData(surface,area,AR))

    return surface
    
end


# Creates beam elements
function create_beam_elements(geometry,sectional,discretization,hinges,attachments)

    # Unpack beam data
    @unpack ℓ,k,R0 = geometry
    @unpack S,I = sectional
    @unpack nElements,normalizedNodalPositions = discretization
    @unpack hingedNodes,hingedNodesDoF = hinges
    @unpack surface = attachments

    # Initialize elements array
    elements = Vector{Element}()

    # Loop to create elements
    for e in 1:nElements
        push!(elements, Element(e-1,ℓ,k,R0,S,I,normalizedNodalPositions,hingedNodes,hingedNodesDoF,surface))
    end

    return elements

end


# Computes the derivatives of the initial generalized displacements with respect to the beam arclength coordinate, x₁
function initial_displacements_derivatives(u0_of_x1,p0_of_x1)

    uprime0_of_x1 = x1 -> ForwardDiff.derivative(u0_of_x1, x1)
    pprime0_of_x1 = x1 -> ForwardDiff.derivative(p0_of_x1, x1)

    return uprime0_of_x1,pprime0_of_x1

end