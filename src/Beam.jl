abstract type BeamElement end


"""
@with_kw mutable struct Beam

Beam composite type

# Fields
- name::String = name of the beam
- connectedNodesThis::Int64 = which node is connected to another beam
- connectedBeams::Union{Nothing, Beam} = to which beam the current beam is connected
- connectedNodesOther::Int64 = to which node of the other beam it is connected
- length::Float64 = total arclength of the beam
- rotationParametrization::String = rotation parametrization for the undeformed geometry definition
- p0::Vector{Float64} = corresponding rotation parameters
- k::Vector{Float64} = initial curvatures of the beam in the undeformed geometry (torsional, flapwise bending, in-plane bending)
- nElements::Int64 = number of elements
- normalizedNodalPositions::Vector{Float64} = nodal positions normalized by the length
- elementRange::Vector{Int64} = current beam's range of elements in the model
- constitutiveRelation::{String} = constitutive relation of the material
- C::Matrix{Float64} = sectional stiffness matrix
- I::Matrix{Float64} = sectional inertia matrix

# Notes
 - Some fields have meaningful default values, others need inputs
"""
@with_kw mutable struct Beam
    # Name and ID
    name::String = ""
    ID::Int64 = 0
    # Connection relative to other beams
    connectedBeams::Union{Nothing,Vector{Beam}} = nothing
    connectedNodesThis::Vector{Int64} = Vector{Int64}()
    connectedNodesOther::Vector{Int64} = Vector{Int64}()
    # Geometry 
    length::Number 
    rotationParametrization::String = "WM"
    p0::Vector{Float64} = zeros(3)
    k::Vector{Float64} = zeros(3)
    R0::Matrix{Float64} = I3
    R_cs::Matrix{Float64} = I3
    # Discretization
    nElements::Int64 
    normalizedNodalPositions::Vector{Float64} = Vector{Float64}()
    # Assembly variables
    elementRange::Vector{Int64} = Vector{Int64}()
    nodeRange::Vector{Int64} = Vector{Int64}()
    r_n::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    # Sectional properties (stiffness and inertia matrices)
    C::Union{Vector{Matrix{Float64}},Vector{Matrix{Int64}}} 
    I::Union{Vector{Matrix{Float64}},Vector{Matrix{Int64}}} = [I6]
    # Elements
    elements::Vector{<:BeamElement} = Vector{Element}()
    # Point inertias
    pointInertias::Vector{PointInertia} = Vector{PointInertia}()
    # Velocity DoFs to update on initial conditions
    velDoFToUpdate::BitVector = trues(6)
    # Hinged nodes and hinged DoF
    hingedNodes::Vector{Int64} = Vector{Int64}()
    hingedNodesDoF::Union{Vector{Vector{Bool}},Vector{BitVector}} = Vector{BitVector}()
    # Initial displacements/rotations/velocities
    u0_of_x1::Union{Vector{<:Number},<:Function,Nothing} = nothing
    p0_of_x1::Union{Vector{<:Number},<:Function,Nothing} = nothing
    udot0_of_x1::Union{Vector{<:Number},<:Function,Nothing} = nothing
    pdot0_of_x1::Union{Vector{<:Number},<:Function,Nothing} = nothing
    # Derivatives of initial generalized displacements with respect to arclength coordinate
    uprime0_of_x1::Function = x1 -> zeros(3)
    pprime0_of_x1::Function = x1 -> zeros(3)
    # Distributed loads
    f_A_of_x1t::Union{Nothing,<:Function} = nothing
    m_A_of_x1t::Union{Nothing,<:Function} = nothing
    f_b_of_x1t::Union{Nothing,<:Function} = nothing
    m_b_of_x1t::Union{Nothing,<:Function} = nothing
    ff_A_of_x1t::Union{Nothing,<:Function} = nothing
    mf_A_of_x1t::Union{Nothing,<:Function} = nothing
    ff_b_of_x1t::Union{Nothing,<:Function} = nothing
    mf_b_of_x1t::Union{Nothing,<:Function} = nothing
    # TFs for non-zero distributed loads
    hasDistributedDeadForcesBasisA::Bool = false
    hasDistributedDeadMomentsBasisA::Bool = false
    hasDistributedDeadForcesBasisb::Bool = false
    hasDistributedDeadMomentsBasisb::Bool = false
    hasDistributedFollowerForcesBasisA::Bool = false
    hasDistributedFollowerMomentsBasisA::Bool = false
    hasDistributedFollowerForcesBasisb::Bool = false
    hasDistributedFollowerMomentsBasisb::Bool = false

    # Beam constructor
    function Beam(name::String,ID::Int64,connectedBeams::Union{Nothing,Vector{Beam}},connectedNodesThis::Vector{Int64},connectedNodesOther::Vector{Int64},length::Number,rotationParametrization::String,p0::Vector{Float64},k::Vector{Float64},R0::Matrix{Float64},R_cs::Matrix{Float64},nElements::Int64,normalizedNodalPositions::Vector{Float64},elementRange::Vector{Int64},nodeRange::Vector{Int64},r_n::Vector{Vector{Float64}},C::Union{Vector{Matrix{Float64}},Vector{Matrix{Int64}}},I::Union{Vector{Matrix{Float64}},Vector{Matrix{Int64}}},elements::Vector{<:BeamElement},pointInertias::Vector{PointInertia},velDoFToUpdate::BitVector,hingedNodes::Vector{Int64},hingedNodesDoF::Union{Vector{Vector{Bool}},Vector{BitVector}},u0_of_x1::Union{Vector{<:Number},<:Function,Nothing},p0_of_x1::Union{Vector{<:Number},<:Function,Nothing},udot0_of_x1::Union{Vector{<:Number},<:Function,Nothing},pdot0_of_x1::Union{Vector{<:Number},<:Function,Nothing},uprime0_of_x1::Function,pprime0_of_x1::Function,f_A_of_x1t::Union{Nothing,<:Function},m_A_of_x1t::Union{Nothing,<:Function},f_b_of_x1t::Union{Nothing,<:Function},m_b_of_x1t::Union{Nothing,<:Function},ff_A_of_x1t::Union{Nothing,<:Function},mf_A_of_x1t::Union{Nothing,<:Function},ff_b_of_x1t::Union{Nothing,<:Function},mf_b_of_x1t::Union{Nothing,<:Function},hasDistributedDeadForcesBasisA::Bool,hasDistributedDeadMomentsBasisA::Bool,hasDistributedDeadForcesBasisb::Bool,hasDistributedDeadMomentsBasisb::Bool,hasDistributedFollowerForcesBasisA::Bool,hasDistributedFollowerMomentsBasisA::Bool,hasDistributedFollowerForcesBasisb::Bool,hasDistributedFollowerMomentsBasisb::Bool) 

        # Initialize the beam
        self = new(name,ID,connectedBeams,connectedNodesThis,connectedNodesOther,length,rotationParametrization,p0,k,R0,R_cs,nElements,normalizedNodalPositions,elementRange,nodeRange,r_n,C,I,elements,pointInertias,velDoFToUpdate,hingedNodes,hingedNodesDoF,u0_of_x1,p0_of_x1,udot0_of_x1,pdot0_of_x1,uprime0_of_x1,pprime0_of_x1,f_A_of_x1t,m_A_of_x1t,f_b_of_x1t,m_b_of_x1t,ff_A_of_x1t,mf_A_of_x1t,ff_b_of_x1t,mf_b_of_x1t,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb)

        # Validate and update the beam 
        update_beam!(self)

        return self
    end
end
export Beam


"""
update_beam!(beam::Beam)

Validates and updates the beam construction

# Arguments
- beam::Beam
"""
function update_beam!(beam::Beam)

    # Validate sectional matrices
    validate_sectional_matrices(beam)

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

    # Get velocity DoFs to update on initial dynamic time step
    get_velocity_dofs_to_update!(beam)

    # Validate rotation parameters and get rotation tensors 
    get_rotation_tensors!(beam)

    # Get derivatives of generalized initial displcements 
    get_initial_displacements_derivatives!(beam)

    # Create beam elements
    create_beam_elements!(beam)

    # Set nodal coordinates
    set_nodal_coordinates!(beam)

    return beam

end
export update_beam!


"""
validate_sectional_matrices(beam::Beam)

Checks that sectional matrices are input as a single one for the whole beam or are input in a per element basis 

# Arguments
- beam::Beam
"""
function validate_sectional_matrices(beam::Beam)

    @unpack nElements,C,I = beam

    # Sectional stiffness matrix
    @assert (length(C)==1 || length(C)==nElements)

    # Sectional inertia matrix
    @assert (length(I)==1 || length(I)==nElements)

end


"""
validate_connected_beams(beam::Beam)

Checks that the beam connections are consistent

# Arguments
- beam::Beam
"""
function validate_connected_beams(beam::Beam)

    @unpack connectedBeams,connectedNodesThis,connectedNodesOther = beam

    if !isnothing(connectedBeams)
        @assert length(connectedBeams) == length(connectedNodesThis) == length(connectedNodesOther) "connectedBeams, connectedNodesThis and connectedNodesOther arrays must have the same length"
    end

end


"""
get_velocity_dofs_to_update!(beam::Beam)

Gets the velocity DOFs (V and Ω) to be updated on the initial dynamic solution

# Arguments
- beam::Beam
"""
function get_velocity_dofs_to_update!(beam::Beam)

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


"""
validate_initial_conditions!(beam::Beam)

Validate initial generalized displacements/velocities defined in basis b, and transform to function if input as constant vector

# Arguments
- beam::Beam
"""
function validate_initial_conditions!(beam::Beam)

    @unpack u0_of_x1,p0_of_x1,udot0_of_x1,pdot0_of_x1 = beam

    # Displacements 
    if u0_of_x1 isa Function
        @assert length(u0_of_x1(0)) == 3
    elseif u0_of_x1 isa Vector{<:Number}
        @assert length(u0_of_x1) == 3
        u0_of_x1_const = deepcopy(u0_of_x1)
        u0_of_x1 = x1 -> u0_of_x1_const
    elseif isnothing(u0_of_x1)
        u0_of_x1 = x1 -> zeros(3)
    end

    # Rotation parameters
    if p0_of_x1 isa Function
        @assert length(p0_of_x1(0)) == 3
    elseif p0_of_x1 isa Vector{<:Number}
        @assert length(p0_of_x1) == 3
        p0_of_x1_const = deepcopy(p0_of_x1)
        p0_of_x1 = x1 -> p0_of_x1_const
    elseif isnothing(p0_of_x1)
        p0_of_x1 = x1 -> zeros(3)    
    end

    # Displacements' rates
    if udot0_of_x1 isa Function
        @assert length(udot0_of_x1(0)) == 3
    elseif udot0_of_x1 isa Vector{<:Number}
        @assert length(udot0_of_x1) == 3
        udot0_of_x1_const = deepcopy(udot0_of_x1)
        udot0_of_x1 = x1 -> udot0_of_x1_const
    elseif isnothing(udot0_of_x1)
        udot0_of_x1 = x1 -> zeros(3)    
    end

    # Rotation parameters' rates
    if pdot0_of_x1 isa Function
        @assert length(pdot0_of_x1(0)) == 3
    elseif pdot0_of_x1 isa Vector{<:Number}
        @assert length(pdot0_of_x1) == 3
        pdot0_of_x1_const = deepcopy(pdot0_of_x1)
        pdot0_of_x1 = x1 -> pdot0_of_x1_const
    elseif isnothing(pdot0_of_x1)
        pdot0_of_x1 = x1 -> zeros(3)    
    end

    @pack! beam = u0_of_x1,p0_of_x1,udot0_of_x1,pdot0_of_x1

end


"""
validate_normalized_nodal_positions!(beam::Beam)

Validates the normalized nodal positions if they were input, or updates them if they were not

# Arguments
- beam::Beam
"""
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


"""
validate_hinged_nodes!(beam::Beam)
 
Validates the hinged nodes and DoFs

# Arguments
- beam::Beam
"""
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


"""
validate_distributed_loads!(beam::Beam)
 
Validates and updates the distributed loads

# Arguments
- beam::Beam
"""
function validate_distributed_loads!(beam::Beam)

    @unpack f_A_of_x1t,m_A_of_x1t,f_b_of_x1t,m_b_of_x1t,ff_A_of_x1t,mf_A_of_x1t,ff_b_of_x1t,mf_b_of_x1t,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb = beam

    if isnothing(f_A_of_x1t)
        f_A_of_x1t = (x1,t) -> zeros(3)
    else
        @assert isa(f_A_of_x1t(0,0),Vector{<:Number})
        @assert length(f_A_of_x1t(0,0)) == 3
        hasDistributedDeadForcesBasisA = true
    end
    if isnothing(m_A_of_x1t)
        m_A_of_x1t = (x1,t) -> zeros(3)
    else
        @assert isa(m_A_of_x1t(0,0),Vector{<:Number})
        @assert length(m_A_of_x1t(0,0)) == 3
        hasDistributedDeadMomentsBasisA = true
    end
    if isnothing(f_b_of_x1t)
        f_b_of_x1t = (x1,t) -> zeros(3)
    else
        @assert isa(f_b_of_x1t(0,0),Vector{<:Number})
        @assert length(f_b_of_x1t(0,0)) == 3
        hasDistributedDeadForcesBasisb = true
    end
    if isnothing(m_b_of_x1t)
        m_b_of_x1t = (x1,t) -> zeros(3)
    else
        @assert isa(m_b_of_x1t(0,0),Vector{<:Number})
        @assert length(m_b_of_x1t(0,0)) == 3
        hasDistributedDeadMomentsBasisb = true
    end
    if isnothing(ff_A_of_x1t)
        ff_A_of_x1t = (x1,t) -> zeros(3)
    else
        @assert isa(ff_A_of_x1t(0,0),Vector{<:Number})
        @assert length(ff_A_of_x1t(0,0)) == 3
        hasDistributedFollowerForcesBasisA = true
    end
    if isnothing(mf_A_of_x1t)
        mf_A_of_x1t = (x1,t) -> zeros(3)
    else
        @assert isa(mf_A_of_x1t(0,0),Vector{<:Number})
        @assert length(mf_A_of_x1t(0,0)) == 3
        hasDistributedFollowerMomentsBasisA = true
    end
    if isnothing(ff_b_of_x1t)
        ff_b_of_x1t = (x1,t) -> zeros(3)
    else
        @assert isa(ff_b_of_x1t(0,0),Vector{<:Number})
        @assert length(ff_b_of_x1t(0,0)) == 3
        hasDistributedFollowerForcesBasisb = true
    end
    if isnothing(mf_b_of_x1t)
        mf_b_of_x1t = (x1,t) -> zeros(3)
    else
        @assert isa(mf_b_of_x1t(0,0),Vector{<:Number})
        @assert length(mf_b_of_x1t(0,0)) == 3
        hasDistributedFollowerMomentsBasisb = true
    end

    @pack! beam = f_A_of_x1t,m_A_of_x1t,f_b_of_x1t,m_b_of_x1t,ff_A_of_x1t,mf_A_of_x1t,ff_b_of_x1t,mf_b_of_x1t,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb
end


"""
get_rotation_tensors!(beam::Beam)

Validates the input rotation parameters and rotation parametrization, and gets the corresponding rotation tensors from basis A to basis b and of the cross-section

# Arguments
- beam::Beam
"""
function get_rotation_tensors!(beam::Beam)

    @unpack p0,rotationParametrization = beam

    # Validate
    @assert in(rotationParametrization,["E321","E313","WM"])
    @assert length(p0) == 3

    # Get rotation tensors
    if rotationParametrization == "E321"
        R0 = rotation_tensor_E321(p0)
        # R_cs = rotation_tensor_E321([0.0,0.0,p0[3]])
    elseif rotationParametrization == "E313"
        R0 = get_rotation_tensor_E313(p0)
        # R_cs = get_rotation_tensor_E313([0.0,p0[3],0.0])
    elseif rotationParametrization == "WM"
        R0,_ = rotation_tensor_WM(p0)
        # R_cs,_ = rotation_tensor_WM([p0[3],0.0,0.0])
    end

    R_cs = I3

    @pack! beam = R0,R_cs
end


"""
get_initial_displacements_derivatives!(beam::Beam)

Gets the derivatives of the initial generalized displacements

# Arguments
- beam::Beam
"""
function get_initial_displacements_derivatives!(beam::Beam)

    @unpack u0_of_x1,p0_of_x1 = beam

    # Get the derivative functions with respect to x1
    uprime0_of_x1 = x1 -> ForwardDiff.derivative(u0_of_x1, x1)
    pprime0_of_x1 = x1 -> ForwardDiff.derivative(p0_of_x1, x1)

    @pack! beam = uprime0_of_x1,pprime0_of_x1

end


"""
create_beam_elements!(beam::Beam)

Initializes the element and node ranges, and creates beam elements

# Arguments
- beam::Beam
"""
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


"""
set_nodal_coordinates!(beam::Beam)

Sets the nodal coordinates of each element into the beam

# Arguments
- beam::Beam
"""
function set_nodal_coordinates!(beam::Beam)

    # Reset nodal coordinates
    r_n = Vector{Vector{Float64}}()

    # Set nodal coordinates of each element into the beam
    for element in beam.elements
        push!(r_n,element.r_n1,element.r_n2)
    end

    # Remove duplicated coordinates
    unique!(r_n)

    @pack! beam = r_n

end


"""
add_point_inertias_to_beam!(;beam::Beam,inputPointInertias::Vector{PointInertia})

Adds loads to the beam

# Arguments
- beam::Beam
- inputPointInertias::Vector{PointInertia}
"""
function add_point_inertias_to_beam!(beam::Beam;inputPointInertias::Vector{PointInertia})

    @unpack nElements,elements,pointInertias = beam

    # Loop point inertias
    for pointInertia in inputPointInertias

        @unpack elementLocalID = pointInertia

        # Check that the beam has the element to which the point inertia was assigned
        @assert elementLocalID <= nElements

        # Add point inertia to beam
        push!(pointInertias,pointInertia)

        # Add point inertia to element (increment its sectional inertia matrix)
        add_point_inertia_to_element!(elements[elementLocalID],pointInertia)
    end

    @pack! beam = pointInertias

end
export add_point_inertias_to_beam!


"""
add_loads_to_beam!(;beam::Beam,loadTypes::Vector{String},loadFuns::Vector{<:Function})

Adds loads to the beam

# Arguments
- beam::Beam
- loadTypes::Vector{String}
- loadFuns::Vector{<:Function}
"""
function add_loads_to_beam!(beam::Beam;loadTypes::Vector{String},loadFuns::Vector{<:Function})
    
    # Loop load types and respective functions
    for (load,fun) in zip(loadTypes,loadFuns)

        # Check inputs
        @assert in(load,["f_A_of_x1t","m_A_of_x1t","f_b_of_x1t","m_b_of_x1t","ff_A_of_x1t","mf_A_of_x1t","ff_b_of_x1t","mf_b_of_x1t"])
        @assert isa(fun(0,0),Vector{<:Number})
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

        # Loop elements
        for element in beam.elements
            update_element_distributed_loads!(element,load,fun)
        end

    end

end
export add_loads_to_beam!


"""
add_initial_displacements_and_velocities_to_beam!(beam::Beam,conditionTypes::Vector{String},conditionFuns::Vector{<:Function})

Adds initial generalized displacements and velocities to the beam

# Arguments
- beam::Beam
- conditionTypes::Vector{String}
- conditionFuns::Vector{<:Function}
"""
function add_initial_displacements_and_velocities_to_beam!(beam::Beam;conditionTypes::Vector{String},conditionFuns::Vector{<:Function})
    
    # Loop initial condition types and respective functions
    for (type,fun) in zip(conditionTypes,conditionFuns)

        # Check inputs
        @assert in(type,["u0_of_x1","p0_of_x1","udot0_of_x1","pdot0_of_x1"])
        @assert isa(fun(0),Vector{<:Number})
        @assert length(fun(0)) == 3

        # Update initial condition on beam
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

end
export add_initial_displacements_and_velocities_to_beam!