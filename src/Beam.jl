abstract type BeamElement end

#

    # Beam composite type

#
@with_kw mutable struct Beam

    # Primary (inputs to beam creation)
    # ---------------------------------
    # Name 
    name::String
    # Geometry 
    length::Real 
    rotationParametrization::String
    p0::Vector{<:Real}
    k::Union{Vector{<:Real},<:Function}
    initialPosition::Vector{<:Real}
    # Discretization
    nElements::Int64 
    normalizedNodalPositions::Vector{Float64}
    # Sectional properties (stiffness and inertia matrices)
    S::Vector{<:Matrix{<:Real}} 
    I::Vector{<:Matrix{<:Real}}
    # Connection relative to other beams
    connectedBeams::Union{Nothing,Vector{Beam}}
    connectedNodesThis::Vector{Int64}
    connectedNodesOther::Vector{Int64}
    # Attached point inertias
    pointInertias::Vector{PointInertia}
    # Hinged nodes and hinged DoF
    hingedNodes::Vector{Int64}
    hingedNodesDoF::Union{Vector{Vector{Bool}},Vector{BitVector}}
    # Initial generalized displacements and velocities
    u0_of_x1::Union{Vector{<:Real},<:Function,Nothing}
    p0_of_x1::Union{Vector{<:Real},<:Function,Nothing}
    udot0_of_x1::Union{Vector{<:Real},<:Function,Nothing}
    pdot0_of_x1::Union{Vector{<:Real},<:Function,Nothing}
    # Distributed loads
    f_A_of_x1t::Union{Nothing,<:Function}
    m_A_of_x1t::Union{Nothing,<:Function}
    f_b_of_x1t::Union{Nothing,<:Function}
    m_b_of_x1t::Union{Nothing,<:Function}
    ff_A_of_x1t::Union{Nothing,<:Function}
    mf_A_of_x1t::Union{Nothing,<:Function}
    ff_b_of_x1t::Union{Nothing,<:Function}
    mf_b_of_x1t::Union{Nothing,<:Function}
    # Attached aerodynamic surface
    aeroSurface::Union{Nothing,AeroSurface}
    # Attached springs
    springs::Vector{Spring}

    # Secondary (outputs from beam creation)
    # --------------------------------------
    # Elements
    elements::Vector{<:BeamElement} = Vector{Element}()
    # Rotation tensor from basis A to basis b
    R0::Matrix{Float64} = I3
    # Assembly variables
    ID::Int64 = 0
    elementRange::Vector{Int64} = Vector{Int64}()
    nodeRange::Vector{Int64} = Vector{Int64}()
    r_n::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    # Velocity DoFs to update on initial conditions
    velDoFToUpdate::BitVector = trues(6)
    # Derivatives of initial generalized displacements with respect to arclength coordinate
    uprime0_of_x1::Function = x1 -> zeros(3)
    pprime0_of_x1::Function = x1 -> zeros(3)
    # TFs for non-zero distributed loads
    hasDistributedDeadForcesBasisA::Bool = false
    hasDistributedDeadMomentsBasisA::Bool = false
    hasDistributedDeadForcesBasisb::Bool = false
    hasDistributedDeadMomentsBasisb::Bool = false
    hasDistributedFollowerForcesBasisA::Bool = false
    hasDistributedFollowerMomentsBasisA::Bool = false
    hasDistributedFollowerForcesBasisb::Bool = false
    hasDistributedFollowerMomentsBasisb::Bool = false

end
export Beam


"""
    create_Beam(; kwargs...)

Beam constructor

# Keyword Arguments
- `name::String`: name of the beam
- `length::Real`: (arc)length of the beam
- `rotationParametrization::String`: type of rotation parametrization to define basis b
- `p0::Vector{<:Real}`: rotation parameters from basis A to basis b
- `k::Union{Vector{<:Real},<:Function}`: undeformed beam's curvatures per unit length
- `initialPosition::Vector{<:Real}`: initial position of the beam's first node relative to the beam's origin (which may be another beam's node)
- `nElements::Int64`: number of elements for discretization
- `normalizedNodalPositions::Vector{Float64}`: normalized nodal positions of beam elements
- `S::Vector{<:Matrix{<:Real}}`: array of sectional stiffness matrices
- `I::Vector{<:Matrix{<:Real}}`: array of sectional inertia matrices
- `connectedBeams::Union{Nothing,Vector{Beam}}`: array of beams to which this beam is connected (a non-recursive property)
- `connectedNodesThis::Vector{Int64}`: nodes of this beam which are connected to other beams' nodes
- `connectedNodesOther::Vector{Int64}`: respective nodes of the other beams
- `pointInertias::Vector{PointInertia}`: attached point inertias
- `hingedNodes::Vector{Int64}`: nodes with a hinge
- `hingedNodesDoF::Union{Vector{Vector{Bool}},Vector{BitVector}}`: respective hinged degrees-of-freedom
- `u0_of_x1::Union{Vector{<:Real},<:Function,Nothing}`: initial displacement (resolved in the undeformed beam basis, b) of the beam as a function of its arclength coordinate (x1)
- `p0_of_x1::Union{Vector{<:Real},<:Function,Nothing}`: initial rotation parameters (resolved in the undeformed beam basis, b) of the beam as a function of its arclength coordinate (x1)
- `udot0_of_x1::Union{Vector{<:Real},<:Function,Nothing}`: initial displacement's rates (resolved in the undeformed beam basis, b) of the beam as a function of its arclength coordinate (x1)
- `pdot0_of_x1::Union{Vector{<:Real},<:Function,Nothing}`: initial rotation parameters' rates (resolved in the undeformed beam basis, b) of the beam as a function of its arclength coordinate (x1)
- `f_A_of_x1t::Union{Nothing,<:Function}`: distributed dead forces initially resolved in basis A, as a function of the beam arclength coordinate (x1) and time (t)
- `m_A_of_x1t::Union{Nothing,<:Function}`: distributed dead moments initially resolved in basis A, as a function of the beam arclength coordinate (x1) and time (t)
- `f_b_of_x1t::Union{Nothing,<:Function}`: distributed dead forces initially resolved in basis b, as a function of the beam arclength coordinate (x1) and time (t)
- `m_b_of_x1t::Union{Nothing,<:Function}`: distributed dead moments initially resolved in basis b, as a function of the beam arclength coordinate (x1) and time (t)
- `ff_A_of_x1t::Union{Nothing,<:Function}`: distributed follower forces initially resolved in basis A, as a function of the beam arclength coordinate (x1) and time (t)
- `mf_A_of_x1t::Union{Nothing,<:Function}`: distributed moments forces initially resolved in basis A, as a function of the beam arclength coordinate (x1) and time (t)
- `ff_b_of_x1t::Union{Nothing,<:Function}`: distributed follower forces initially resolved in basis b, as a function of the beam arclength coordinate (x1) and time (t)
- `mf_b_of_x1t::Union{Nothing,<:Function}`: distributed follower moments initially resolved in basis b, as a function of the beam arclength coordinate (x1) and time (t)
- `aeroSurface::Union{Nothing,AeroSurface}`: attached aerodynamic surface
- `springs::Vector{Spring}`: array of attached springs
"""
function create_Beam(; name::String="",length::Real,rotationParametrization::String="E321",p0::Vector{<:Real}=zeros(3),k::Union{Vector{<:Real},<:Function}=zeros(3),initialPosition::Vector{<:Real}=zeros(3),nElements::Int64,normalizedNodalPositions::Vector{Float64}=Vector{Float64}(),S::Vector{<:Matrix{<:Real}},I::Vector{<:Matrix{<:Real}}=[I6],connectedBeams::Union{Nothing,Vector{Beam}}=nothing,connectedNodesThis::Vector{Int64}=Vector{Int64}(),connectedNodesOther::Vector{Int64}=Vector{Int64}(),pointInertias::Vector{PointInertia}=Vector{PointInertia}(),hingedNodes::Vector{Int64}=Vector{Int64}(),hingedNodesDoF::Union{Vector{Vector{Bool}},Vector{BitVector}}=Vector{BitVector}(),u0_of_x1::Union{Vector{<:Real},<:Function,Nothing}=nothing,p0_of_x1::Union{Vector{<:Real},<:Function,Nothing}=nothing,udot0_of_x1::Union{Vector{<:Real},<:Function,Nothing}=nothing,pdot0_of_x1::Union{Vector{<:Real},<:Function,Nothing}=nothing,f_A_of_x1t::Union{Nothing,<:Function}=nothing,m_A_of_x1t::Union{Nothing,<:Function}=nothing,f_b_of_x1t::Union{Nothing,<:Function}=nothing,m_b_of_x1t::Union{Nothing,<:Function}=nothing,ff_A_of_x1t::Union{Nothing,<:Function}=nothing,mf_A_of_x1t::Union{Nothing,<:Function}=nothing,ff_b_of_x1t::Union{Nothing,<:Function}=nothing,mf_b_of_x1t::Union{Nothing,<:Function}=nothing,aeroSurface::Union{Nothing,AeroSurface}=nothing,springs::Vector{Spring}=Vector{Spring}())

    # Initialize the beam
    self = Beam(name=name,length=length,rotationParametrization=rotationParametrization,p0=p0,k=k,initialPosition=initialPosition,nElements=nElements,normalizedNodalPositions=normalizedNodalPositions,S=S,I=I,connectedBeams=connectedBeams,connectedNodesThis=connectedNodesThis,connectedNodesOther=connectedNodesOther,pointInertias=pointInertias,hingedNodes=hingedNodes,hingedNodesDoF=hingedNodesDoF,u0_of_x1=u0_of_x1,p0_of_x1=p0_of_x1,udot0_of_x1=udot0_of_x1,pdot0_of_x1=pdot0_of_x1,f_A_of_x1t=f_A_of_x1t,m_A_of_x1t=m_A_of_x1t,f_b_of_x1t=f_b_of_x1t,m_b_of_x1t=m_b_of_x1t,ff_A_of_x1t=ff_A_of_x1t,mf_A_of_x1t=mf_A_of_x1t,ff_b_of_x1t=ff_b_of_x1t,mf_b_of_x1t=mf_b_of_x1t,aeroSurface=aeroSurface,springs=springs)

    # Validate and update the beam 
    update_beam!(self)

    return self
end
export create_Beam


"""
    update_beam!(beam::Beam)

Validates and updates the beam construction

# Arguments
- `beam::Beam`
"""
function update_beam!(beam::Beam)

    # Validate beam
    validate_beam!(beam)

    # Get velocity DoFs to update on initial dynamic time step
    velocity_dofs_to_update!(beam)

    # Get rotation tensor from basis A to basis b 
    get_rotation_tensor!(beam)

    # Get derivatives of generalized initial displacements with respect to arclength coordinate
    initial_displacements_derivatives!(beam)

    # Update attached aerodynamic surface properties
    update_aero_surface!(beam)

    # Create beam elements
    create_beam_elements!(beam)

    # Set nodal coordinates
    set_nodal_coordinates!(beam)

    return beam

end
export update_beam!


# Validates the beam inputs
function validate_beam!(beam::Beam)

    # Validate initial position of the first node of the beam
    @assert length(beam.initialPosition) == 3

    # Validate sectional matrices
    validate_sectional_matrices(beam)

    # Validate rotation parametrization
    validate_rotation_parametrization(beam)

    # Validate beam connections
    validate_connected_beams(beam)
    
    # Validate initial generalized displacements/velocities 
    validate_initial_conditions!(beam)

    # Validate and update normalized nodal positions
    validate_normalized_nodal_positions!(beam)

    # Validate hinged nodes and hinged DoF
    validate_hinged_nodes!(beam)

    # Validate and update distributed loads inputs
    validate_distributed_loads!(beam)

end


# Checks that sectional matrices are input as a single one for the whole beam or are input in a per element basis 
function validate_sectional_matrices(beam::Beam)

    @unpack nElements,S,I = beam

    # Sectional stiffness matrix
    @assert (length(S)==1 || length(S)==nElements) "input either one stiffness matrix for the entire beam, or one for each element"
    for Si in S
        @assert size(Si) == (6,6) "stiffness matrices must be of size (6,6)"
        @assert all([Si[j,j] > 0 for j in 1:6]) "diagonal elements of stiffness matrices must be positive"
    end

    # Sectional inertia matrix
    @assert (length(I)==1 || length(I)==nElements) "input either one inertia matrix for the entire beam, or one for each element"
    for Ii in I
        @assert size(Ii) == (6,6) "inertia matrices must be of size (6,6)"
        @assert all([Ii[j,j] >= 0 for j in 1:6]) "diagonal elements of inertia matrices must be greater than or equal to zero"
    end

end


# Validates the input rotation parameters and rotation parametrization
function validate_rotation_parametrization(beam::Beam)

    @unpack p0,rotationParametrization = beam

    @assert rotationParametrization in ["E321","E213","E231","E313"]
    @assert length(p0) == 3

end


# Checks that the beam connections are consistent
function validate_connected_beams(beam::Beam)

    @unpack connectedBeams,connectedNodesThis,connectedNodesOther = beam

    if !isnothing(connectedBeams)
        @assert length(connectedBeams) == length(connectedNodesThis) == length(connectedNodesOther) "connectedBeams, connectedNodesThis and connectedNodesOther arrays must have the same length"
    end

end


# Validate initial generalized displacements/velocities defined in basis b, and transform to function if input as constant vector
function validate_initial_conditions!(beam::Beam)

    @unpack u0_of_x1,p0_of_x1,udot0_of_x1,pdot0_of_x1 = beam

    # Displacements 
    if u0_of_x1 isa Function
        @assert length(u0_of_x1(0)) == 3
    elseif u0_of_x1 isa Vector{<:Real}
        @assert length(u0_of_x1) == 3
        u0_of_x1_const = deepcopy(u0_of_x1)
        u0_of_x1 = x1 -> u0_of_x1_const
    elseif isnothing(u0_of_x1)
        u0_of_x1 = x1 -> zeros(3)
    end

    # Rotation parameters
    if p0_of_x1 isa Function
        @assert length(p0_of_x1(0)) == 3
    elseif p0_of_x1 isa Vector{<:Real}
        @assert length(p0_of_x1) == 3
        p0_of_x1_const = deepcopy(p0_of_x1)
        p0_of_x1 = x1 -> p0_of_x1_const
    elseif isnothing(p0_of_x1)
        p0_of_x1 = x1 -> zeros(3)    
    end

    # Displacements' rates
    if udot0_of_x1 isa Function
        @assert length(udot0_of_x1(0)) == 3
    elseif udot0_of_x1 isa Vector{<:Real}
        @assert length(udot0_of_x1) == 3
        udot0_of_x1_const = deepcopy(udot0_of_x1)
        udot0_of_x1 = x1 -> udot0_of_x1_const
    elseif isnothing(udot0_of_x1)
        udot0_of_x1 = x1 -> zeros(3)    
    end

    # Rotation parameters' rates
    if pdot0_of_x1 isa Function
        @assert length(pdot0_of_x1(0)) == 3
    elseif pdot0_of_x1 isa Vector{<:Real}
        @assert length(pdot0_of_x1) == 3
        pdot0_of_x1_const = deepcopy(pdot0_of_x1)
        pdot0_of_x1 = x1 -> pdot0_of_x1_const
    elseif isnothing(pdot0_of_x1)
        pdot0_of_x1 = x1 -> zeros(3)    
    end

    @pack! beam = u0_of_x1,p0_of_x1,udot0_of_x1,pdot0_of_x1

end


# Validates the normalized nodal positions if they were input, or updates them if they were not
function validate_normalized_nodal_positions!(beam::Beam)

    @unpack normalizedNodalPositions,nElements = beam

    if isempty(normalizedNodalPositions)
        normalizedNodalPositions = LinRange(0.0,1.0,nElements+1)
    else
        @assert first(normalizedNodalPositions) == 0.0
        @assert last(normalizedNodalPositions) == 1.0
        @assert length(normalizedNodalPositions) == nElements+1
        @assert issorted(normalizedNodalPositions)
    end

    @pack! beam = normalizedNodalPositions
end


# Validates the hinged nodes and DoFs
function validate_hinged_nodes!(beam::Beam)

    @unpack hingedNodes,hingedNodesDoF,nElements = beam

    lastHingedNode = isempty(hingedNodes) ? 0 : maximum(hingedNodes)
    @assert lastHingedNode <= nElements+1 "last hingedNode greater than number of elements + 1"
    @assert !any(hingedNodes[i]+1==hingedNodes[i+1] for i in 1:length(hingedNodes)-1) "cannot have two consecutive nodes hinged"
    for (i,hingedNodeDoF) in enumerate(hingedNodesDoF)
        @assert length(hingedNodeDoF) == 3 "each vector of hingedNodesDoF must have a length of 3"
        hingedNodesDoF[i] = BitVector(hingedNodeDoF)
    end

    @pack! beam = hingedNodesDoF

end


# Validates and updates the distributed loads
function validate_distributed_loads!(beam::Beam)

    @unpack f_A_of_x1t,m_A_of_x1t,f_b_of_x1t,m_b_of_x1t,ff_A_of_x1t,mf_A_of_x1t,ff_b_of_x1t,mf_b_of_x1t,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb = beam

    if !isnothing(f_A_of_x1t)
        @assert isa(f_A_of_x1t(0,0),Vector{<:Real})
        @assert length(f_A_of_x1t(0,0)) == 3
        hasDistributedDeadForcesBasisA = true
    end
    if !isnothing(m_A_of_x1t)
        @assert isa(m_A_of_x1t(0,0),Vector{<:Real})
        @assert length(m_A_of_x1t(0,0)) == 3
        hasDistributedDeadMomentsBasisA = true
    end
    if !isnothing(f_b_of_x1t)
        @assert isa(f_b_of_x1t(0,0),Vector{<:Real})
        @assert length(f_b_of_x1t(0,0)) == 3
        hasDistributedDeadForcesBasisb = true
    end
    if !isnothing(m_b_of_x1t)
        @assert isa(m_b_of_x1t(0,0),Vector{<:Real})
        @assert length(m_b_of_x1t(0,0)) == 3
        hasDistributedDeadMomentsBasisb = true
    end
    if !isnothing(ff_A_of_x1t)
        @assert isa(ff_A_of_x1t(0,0),Vector{<:Real})
        @assert length(ff_A_of_x1t(0,0)) == 3
        hasDistributedFollowerForcesBasisA = true
    end
    if !isnothing(mf_A_of_x1t)
        @assert isa(mf_A_of_x1t(0,0),Vector{<:Real})
        @assert length(mf_A_of_x1t(0,0)) == 3
        hasDistributedFollowerMomentsBasisA = true
    end
    if !isnothing(ff_b_of_x1t)
        @assert isa(ff_b_of_x1t(0,0),Vector{<:Real})
        @assert length(ff_b_of_x1t(0,0)) == 3
        hasDistributedFollowerForcesBasisb = true
    end
    if !isnothing(mf_b_of_x1t)
        @assert isa(mf_b_of_x1t(0,0),Vector{<:Real})
        @assert length(mf_b_of_x1t(0,0)) == 3
        hasDistributedFollowerMomentsBasisb = true
    end

    @pack! beam = f_A_of_x1t,m_A_of_x1t,f_b_of_x1t,m_b_of_x1t,ff_A_of_x1t,mf_A_of_x1t,ff_b_of_x1t,mf_b_of_x1t,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb

end


# Gets the velocity DOFs (V and Ω) to be updated on the initial dynamic solution
function velocity_dofs_to_update!(beam::Beam)

    @unpack u0_of_x1,p0_of_x1,udot0_of_x1,pdot0_of_x1,length = beam

    # Length discretization
    N = 101
    x1 = LinRange(0,length,N)

    # Initialize
    velDoFToUpdate = trues(6)

    # Loop directions of initial linear velocities
    for i=1:3
        # Get values in current direction
        udot0_of_x1_i = [udot0_of_x1.(x1)[j][i] for j in 1:N]
        # Set flag to update velocity in that direction only if it was input as zero (not input by user)
        velDoFToUpdate[i] = any(!iszero,vcat(udot0_of_x1_i...)) ? false : true
    end

    # Loop directions of initial angular velocities
    for i=4:6
        # Get values in current direction
        pdot0_of_x1_i = [pdot0_of_x1.(x1)[j][i-3] for j in 1:N]
        # Set flag to update angular velocity in that direction only if it was input as zero (not input by user)
        velDoFToUpdate[i] = any(!iszero,vcat(pdot0_of_x1_i...)) ? false : true
    end

    @pack! beam = velDoFToUpdate

end


# Validates the input rotation parameters and rotation parametrization, and gets the corresponding rotation tensor from basis A to basis b 
function get_rotation_tensor!(beam::Beam)

    @unpack p0,rotationParametrization = beam

    # Get rotation tensor
    if rotationParametrization == "E321"
        R0 = rotation_tensor_E321(p0)
    elseif rotationParametrization == "E213"
        R0 = rotation_tensor_E213(p0)
    elseif rotationParametrization == "E231"
        R0 = rotation_tensor_E231(p0)    
    elseif rotationParametrization == "E313"
        R0 = get_rotation_tensor_E313(p0)
    end

    @pack! beam = R0
end


# Computes the derivatives of the initial generalized displacements with respect to the beam arclength coordinate, x₁
function initial_displacements_derivatives!(beam::Beam)

    @unpack u0_of_x1,p0_of_x1 = beam

    uprime0_of_x1 = x1 -> ForwardDiff.derivative(u0_of_x1, x1)
    pprime0_of_x1 = x1 -> ForwardDiff.derivative(p0_of_x1, x1)

    @pack! beam = uprime0_of_x1,pprime0_of_x1

end


# Updates the beam-dependent properties of the attached aerodynamic surface (if any)
function update_aero_surface!(beam::Beam)

    @unpack aeroSurface = beam

    # Skip if there are no attached surface
    if isnothing(aeroSurface)
        return
    end

    # Unpack geometric variables
    @unpack length = beam
    @unpack c,hasSymmetricCounterpart = aeroSurface

    # Compute surface planform area
    if c isa Real
        area = c * length
    else
        area,_ = quadgk(c, 0, length)
    end

    # Aspect ratio
    AR = hasSymmetricCounterpart ? 2*length^2/area : length^2/area

    @pack! beam.aeroSurface = area,AR
end


# Initializes the element and node ranges, and creates beam elements
function create_beam_elements!(beam::Beam)

    @unpack nElements = beam

    # Initialize the element and node range (tied to the model)
    elementRange = 1:nElements
    nodeRange = 1:nElements+1

    # Reset elements
    elements = Vector{Element}()

    @pack! beam = elements,elementRange,nodeRange

    # Create beam elements and assign them to the beam
    for _ in elementRange
        Element(beam)
    end

end


# Sets the nodal coordinates of each element into the beam
function set_nodal_coordinates!(beam::Beam)

    # Reset nodal coordinates
    r_n = Vector{Vector{Float64}}()

    # Set nodal coordinates of each element into the beam (excluding the initial position)
    for element in beam.elements
        push!(r_n,element.r_n1,element.r_n2)
    end

    # Remove duplicated coordinates
    unique!(r_n)

    @pack! beam = r_n

end


"""
    add_point_inertias_to_beam!(beam::Beam; inertias::Vector{PointInertia})

Adds point inertias to the beam

# Arguments
- `beam::Beam`

# Keyword arguments
- `inertias::Vector{PointInertia}`
"""
function add_point_inertias_to_beam!(beam::Beam; inertias::Vector{PointInertia})

    @unpack nElements,elements,pointInertias = beam

    # Loop point inertias
    for pointInertia in inertias

        @unpack elementID = pointInertia

        # Check that the beam has the element to which the point inertia was assigned
        @assert elementID <= nElements

        # Add point inertia to beam
        push!(pointInertias,pointInertia)

        # Add point inertia to element (increment its sectional inertia matrix)
        add_point_inertia_to_element!(elements[elementID],pointInertia)
    end

    @pack! beam = pointInertias

end
export add_point_inertias_to_beam!


"""
    add_loads_to_beam!(beam::Beam; loadTypes::Vector{String},loadFuns::Vector{<:Function})

Adds loads to the beam

# Arguments
- `beam::Beam`

# Keyword arguments
- `loadTypes::Vector{String}`
- `loadFuns::Vector{<:Function}`
"""
function add_loads_to_beam!(beam::Beam; loadTypes::Vector{String},loadFuns::Vector{<:Function})
    
    # Loop load types and respective functions
    for (load,fun) in zip(loadTypes,loadFuns)

        # Check inputs
        @assert in(load,["f_A_of_x1t","m_A_of_x1t","f_b_of_x1t","m_b_of_x1t","ff_A_of_x1t","mf_A_of_x1t","ff_b_of_x1t","mf_b_of_x1t"])
        @assert isa(fun(0,0),Vector{<:Real})
        @assert length(fun(0,0)) == 3

        # Update load on beam
        if load == "f_A_of_x1t"
            beam.f_A_of_x1t = fun
            beam.hasDistributedDeadForcesBasisA = true
        elseif load == "m_A_of_x1t"
            beam.m_A_of_x1t = fun
            beam.hasDistributedDeadMomentsBasisA = true
        elseif load == "f_b_of_x1t"
            beam.f_b_of_x1t = fun
            beam.hasDistributedDeadForcesBasisb = true
        elseif load == "m_b_of_x1t"
            beam.m_b_of_x1t = fun
            beam.hasDistributedDeadMomentsBasisb = true
        elseif load == "ff_A_of_x1t"
            beam.ff_A_of_x1t = fun
            beam.hasDistributedFollowerForcesBasisA = true
        elseif load == "mf_A_of_x1t"
            beam.mf_A_of_x1t = fun
            beam.hasDistributedFollowerMomentsBasisA = true
        elseif load == "ff_b_of_x1t"
            beam.ff_b_of_x1t = fun
            beam.hasDistributedFollowerForcesBasisb = true
        elseif load == "mf_b_of_x1t"
            beam.mf_b_of_x1t = fun
            beam.hasDistributedFollowerMomentsBasisb = true
        end

    end

    # Check validity
    validate_distributed_loads!(beam)

    # Loop load types and respective functions
    for (load,fun) in zip(loadTypes,loadFuns)
        # Loop elements
        for element in beam.elements
            update_element_distributed_loads!(element,load,fun)
        end
    end

end
export add_loads_to_beam!


"""
    add_initial_displacements_and_velocities_to_beam!(beam::Beam; conditionTypes::Vector{String},conditionFuns::Vector{<:Function})

Adds initial generalized displacements and velocities to the beam

# Arguments
- `beam::Beam`

# Keyword arguments
- `conditionTypes::Vector{String}`
- `conditionFuns::Vector{<:Function}`
"""
function add_initial_displacements_and_velocities_to_beam!(beam::Beam;conditionTypes::Vector{String},conditionFuns::Vector{<:Function})
    
    # Loop initial condition types and respective functions
    for (type,fun) in zip(conditionTypes,conditionFuns)

        # Check inputs
        @assert in(type,["u0_of_x1","p0_of_x1","udot0_of_x1","pdot0_of_x1"])
        @assert isa(fun(0),Vector{<:Real})
        @assert length(fun(0)) == 3

        # Set initial condition on beam
        if type == "u0_of_x1"
            beam.u0_of_x1 = fun
        elseif type == "p0_of_x1"
            beam.p0_of_x1 = fun
        elseif type == "udot0_of_x1"
            beam.udot0_of_x1 = fun
        elseif type == "pdot0_of_x1"
            beam.pdot0_of_x1 = fun
        end

        # The update of elemental states occurs on the model creation/update (because V and Ω depend on the model's v_A and ω_A)
    end

    # Check validity  
    validate_initial_conditions!(beam)

    # Get derivatives of generalized initial displacements 
    initial_displacements_derivatives!(beam)

end
export add_initial_displacements_and_velocities_to_beam!


"""
    add_springs_to_beam!(; beam::Beam,springs::Vector{Spring})

Adds simply-attached springs to a beam

# Keyword arguments
- `beam::Beam`
- `springs::Vector{Spring}`
"""
function add_springs_to_beam!(; beam::Beam,springs::Vector{Spring})

    # Loop springs
    for spring in springs
        # Check that the spring is not doubly attached
        @assert !spring.hasDoubleAttachment
        # Check that the beam has the element to which the spring was assigned 
        @assert spring.elementsIDs[1] <= beam.nElements
        # Add spring to beam
        push!(beam.springs,spring)
    end

end
export add_springs_to_beam!


"""
    add_spring_to_beams!(; beams::Vector{Beam},spring::Spring)

Adds a doubly-attached spring to the beams

# Keyword arguments
- `beams::Vector{Beam}`
- `spring::Spring`
"""
function add_spring_to_beams!(; beams::Vector{Beam},spring::Spring)

    # Check that the spring has double attachment and that two beams were input
    @assert spring.hasDoubleAttachment
    @assert length(beams) == 2

    # Loop beams
    for (i,beam) in enumerate(beams)
        # Check that the beam has the element to which the spring was assigned
        @assert spring.elementsIDs[i] <= beam.nElements
        # Add spring to beam
        push!(beam.springs,spring)
    end

end
export add_spring_to_beams!


"""
    remove_all_springs_from_beam!(; beams::Vector{Beam})

Removes all springs attached to the beams

# Keyword arguments
- `beams::Vector{Beam}`
"""
function remove_all_springs_from_beams!(; beams::Vector{Beam})

    for beam in beams
        beam.springs = Vector{Spring}()
    end

end
export remove_all_springs_from_beams!