## BOUNDARY CONDITIONS COMPOSITE TYPES
# ------------------------------------------------------------------------------
# Boundary condition definition composite type
# Fields
# - `name`: name of the BC
# - `beam`: beam at which the BC is applied
# - `node`: local node of the beam at which the BC is applied
# - `types`: types of BCs applied to the node (generalized forces and displacements)
# - `values`: corresponding values of the applied BCs
# - `toBeTrimmed`: boolean vector on whether the BC is to be trimmed
#
struct BCDefinition

    name::String
    beam::Beam
    node::Int
    types::Vector{Symbol}
    values::Vector{<:Function}
    toBeTrimmed::BitVector

end


# Boundary condition deduced data composite type
# Fields
# - `R0_n`: nodal rotation tensor from basis A to basis b
# - `isLoad`: boolean vector for BC type being a load
# - `isFollower`: boolean vector for BC type being a follower load
# - `isTrim`: boolean vector for BC type being a trim load
#
struct BCDeducedData

    R0_n::SMatrix{3,3,Float64,9}
    isLoad::SVector{6,Bool}
    isFollower::SVector{6,Bool}
    isTrim::SVector{6,Bool}

end


# Boundary condition load functions composite type
# Fields
# - `displacementsOnAofTime`: functions of time for displacements on basis A
# - `deadLoadsOnAofTime`: functions of time for dead loads on basis A
# - `followerLoadsOnAofTime`: functions of time for follower loads on basis A
# - `displacementsOnbofTime`: functions of time for displacements on basis b
# - `deadLoadsOnbofTime`: functions of time for dead loads on basis b
# - `followerLoadsOnbofTime`: functions of time for follower loads on basis b
#
struct BCLoadFuns

    displacementsOnAofTime::Vector{Function}
    deadLoadsOnAofTime::Vector{Function}
    followerLoadsOnAofTime::Vector{Function}
    displacementsOnbofTime::Vector{Function}
    deadLoadsOnbofTime::Vector{Function}
    followerLoadsOnbofTime::Vector{Function}

end


# Boundary condition state composite type
# Fields
# - `displacementsOnA`: current generalized displacements on basis A
# - `deadLoadsOnA`: current generalized dead loads on basis A
# - `followerLoadsOnA`: current generalized follower loads on basis A
# - `displacementsOnb`: current generalized displacements on basis b
# - `deadLoadsOnb`: current generalized dead loads on basis b
# - `followerLoadsOnb`: current generalized follower loads on basis b
# - `currentValue`: current generalized displacements/loads vector (displacements and loads are mutually exclusive)
#
struct BCState

    displacementsOnA::MVector{6,Float64}
    deadLoadsOnA::MVector{6,Float64}
    followerLoadsOnA::MVector{6,Float64}
    displacementsOnb::MVector{6,Float64}
    deadLoadsOnb::MVector{6,Float64}
    followerLoadsOnb::MVector{6,Float64}
    currentValue::MVector{6,Float64}

end
BCState() = BCState(
    ntuple(_ -> zeros(MVector{6,Float64}), 7)...
    )


# Boundary condition state cache composite type
# Fields
# - `F̂`: current prescribed force vector resolved in basis A
# - `M̂`: current prescribed moment vector resolved in basis A
# - `Ftrim`: current trim force vector resolved in basis A
# - `Mtrim`: current trim moment vector resolved in basis A
# - `Ff`: current follower force vector resolved in basis A
# - `Mf`: current follower moment vector resolved in basis A
#    
struct BCStateCache

    F̂::MVector{3,Float64}
    M̂::MVector{3,Float64}
    
    Ftrim::MVector{3,Float64}
    Mtrim::MVector{3,Float64}

    Ff::MVector{3,Float64}
    Mf::MVector{3,Float64}

end
BCStateCache() = BCStateCache(
    ntuple(_ -> zeros(MVector{3,Float64}), 2)...,
    ntuple(_ -> zeros(MVector{3,Float64}), 2)...,
    ntuple(_ -> zeros(MVector{3,Float64}), 2)...,
    )  


# Boundary condition history composite type
# Fields
# - `displacementsOnA`: generalized displacements vector over time on basis A
# - `deadLoadsOnA`: generalized dead loads vector over time on basis A
# - `followerLoadsOnA`: generalized follower loads vector over time on basis A
# - `displacementsOnb`: generalized displacements vector over time on basis b
# - `deadLoadsOnb`: generalized dead loads vector over time on basis b
# - `followerLoadsOnb`: generalized follower loads vector over time on basis b
# - `currentValue`: current generalized displacements/loads vector over time
# - `Fmax`: maximum absolute force over time
# - `Mmax`: maximum absolute moment over time
#
struct BCHistory

    displacementsOnA::Vector{MVector{6,Float64}}
    deadLoadsOnA::Vector{MVector{6,Float64}}
    followerLoadsOnA::Vector{MVector{6,Float64}}
    displacementsOnb::Vector{MVector{6,Float64}}
    deadLoadsOnb::Vector{MVector{6,Float64}}
    followerLoadsOnb::Vector{MVector{6,Float64}}
    currentValue::Vector{MVector{6,Float64}}

    Fmax::Float64
    Mmax::Float64
    
end
BCHistory(state::BCState) = BCHistory(
        [state.displacementsOnA],
        [state.deadLoadsOnA],
        [state.followerLoadsOnA],
        [state.displacementsOnb],
        [state.deadLoadsOnb],
        [state.followerLoadsOnb],
        [state.currentValue],
        0.0, 0.0,
    )


# Boundary condition composite type
# Fields
# - `definition`: definition of the BC
# - `deducedData`: deduced data of the BC
# - `loadFuns`: functions of time for the BC
# - `state`: state of the BC
# - `history`: history of the BC over time
#
@with_kw mutable struct BC

    # Definition field (inputs for BC creation)
    definition::BCDefinition

    # Deduced data field (deduced from the definition)
    deducedData::BCDeducedData

    # Loads/displacements functions
    loadFuns::BCLoadFuns

    # State and state cache
    state::BCState
    stateCache::BCStateCache = BCStateCache()

    # History
    history::BCHistory = BCHistory(state)

end
export BC

# Flag for default BCs (no BCs applied)
const _emptyBCflags = [(SVector{6,Bool}(true,true,true,true,true,true), zeros(SVector{6,Bool}), zeros(SVector{6,Bool}))]


## PUBLIC API
# -----------------------------------------------------------------------------
"""
    create_BC(; kwargs...)

BC constructor

# Keyword arguments
- `name::String`: name of the BC
- `beam::Beam`: beam at which the BC is applied
- `node::Int`: node of the beam at which the BC is applied
- `conditions::Dict{Symbol,<:Any}`: BC definition where each key is a BC type symbol and each value is either a scalar constant or a function of time returning a scalar
- `toBeTrimmed::AbstractVector{Bool}`: flag vector on whether the BC is to be trimmed
"""
function create_BC(; 
    name::String="",
    beam::Beam,
    node::Int,
    conditions::Dict{Symbol,<:Any},
    toBeTrimmed::AbstractVector{Bool}=falses(length(conditions)),
    )

    # Validate inputs
    types,values,toBeTrimmed = validate_BC(beam,node,conditions,toBeTrimmed)

    # Initialize BC data
    deducedData,loadFuns,state = initialize_BC_data(beam,node,types,values,toBeTrimmed)

    # Initialize definition field
    definition = BCDefinition(name,beam,node,types,values,toBeTrimmed)

    return BC(
        definition=definition,
        deducedData=deducedData,
        loadFuns=loadFuns,
        state=state
        )

end
export create_BC


"""
    reset_BC_beam!(bc::BC, beam::Beam)

Updates the beam of an existing BC and recomputes beam-dependent BC data

# Arguments
- `bc::BC`: BC to update
- `beam::Beam`: new beam associated with the BC
"""
function reset_BC_beam!(bc::BC,beam::Beam)

    # Unpack current BC definition
    @unpack name,node,types,values,toBeTrimmed = bc.definition

    # Validate with respect to the new beam
    @assert 1 <= node <= beam.discretization.nElements+1 "the BC'ed beam does not contain the BC'ed node"

    # Reset definition beam and recompute deduced data/state for the new beam
    definition = BCDefinition(name,beam,node,types,values,toBeTrimmed)
    deducedData,loadFuns,state = initialize_BC_data(beam,node,types,values,toBeTrimmed)
    history = BCHistory(state)

    # Update BC fields
    setfield!(bc, :definition, definition)
    setfield!(bc, :deducedData, deducedData)
    setfield!(bc, :loadFuns, loadFuns)
    setfield!(bc, :state, state)
    setfield!(bc, :history, history)

end
export reset_BC_beam!


## PRIVATE API
# -----------------------------------------------------------------------------
# Validate input BC data
function validate_BC(beam,node,conditions,toBeTrimmed)

    # Allowed BC types
    allowedTypes = [:u1A,:u2A,:u3A,:p1A,:p2A,:p3A,:u1b,:u2b,:u3b,:p1b,:p2b,:p3b,:F1A,:F2A,:F3A,:M1A,:M2A,:M3A,:F1b,:F2b,:F3b,:M1b,:M2b,:M3b,:Ff1A,:Ff2A,:Ff3A,:Mf1A,:Mf2A,:Mf3A,:Ff1b,:Ff2b,:Ff3b,:Mf1b,:Mf2b,:Mf3b]

    # Initialize BC types and values
    types = Symbol[]
    values = Any[]

    # Set types and values
    for type in allowedTypes
        # Skip if current type is not specified in the conditions dictionary
        haskey(conditions, type) || continue
        # Check if the value is a scalar or a function of time
        valueOrFun = conditions[type]
        if valueOrFun isa Function
            @assert applicable(valueOrFun, 0) "function for $type must be callable with one argument (t)"
            @assert valueOrFun(0) isa Real "$type must return a scalar when called with one argument (t)"
            push!(types, type)
            push!(values, valueOrFun)
        else
            @assert valueOrFun isa Real "$type must be either a scalar or a function of time returning a scalar"
            push!(types, type)
            push!(values, Float64(valueOrFun))
        end
    end

    # Set values as functions of time
    values = map(values) do v
        v isa Function ? 
        v :
        (t -> v)
    end

    # Set flags
    containsDisplacementA = any([type in [:u1A,:u2A,:u3A,:p1A,:p2A,:p3A] for type in types])
    containsDisplacementb = any([type in [:u1b,:u2b,:u3b,:p1b,:p2b,:p3b] for type in types])
    containsDeadLoadA = any([type in [:F1A,:F2A,:F3A,:M1A,:M2A,:M3A] for type in types])
    containsDeadLoadb = any([type in [:F1b,:F2b,:F3b,:M1b,:M2b,:M3b] for type in types])
    containsFollowerLoadA = any([type in [:Ff1A,:Ff2A,:Ff3A,:Mf1A,:Mf2A,:Mf3A] for type in types])
    containsFollowerLoadb = any([type in [:Ff1b,:Ff2b,:Ff3b,:Mf1b,:Mf2b,:Mf3b] for type in types])

    # Check inputs
    @assert all(k -> k in allowedTypes, keys(conditions)) "invalid BC dictionary key, must be one of: $(allowedTypes)"
    @assert types == unique(types) "there are repeated BC types (keys)"
    @assert 1 <= node <= beam.discretization.nElements+1 "the BC'ed beam does not contain the BC'ed node"
    @assert length(types) == length(toBeTrimmed) "set a toBeTrimmed value for each BC type (key)"
    for (i,type) in enumerate(types)
        @assert !(toBeTrimmed[i] && type in [:u1A,:u2A,:u3A,:p1A,:p2A,:p3A,:u1b,:u2b,:u3b,:p1b,:p2b,:p3b]) "Trim variables must be loads"
    end
    containsDisplacementA && containsDisplacementb && error("Specify displacements and rotations either solely on basis A or solely on basis b")
    containsDeadLoadA && containsDeadLoadb && error("Please input dead loads on bases A and b on separate BCs")
    containsFollowerLoadA && containsFollowerLoadb && error("Please input follower loads on bases A and b on separate BCs")

    # Set correct types
    toBeTrimmed = BitVector(toBeTrimmed)

    return types,values,toBeTrimmed

end


# Initializes the BC data
function initialize_BC_data(beam,node,types,values,toBeTrimmed)

    # Get nodal rotation tensor from basis A to basis b
    R0_n = node == length(beam.elements)+1 ?
           beam.elements[node-1].geometry.R0_n2 :
           beam.elements[node].geometry.R0_n1

    # Initialize flags for being loads, being follower loads and being trim loads
    isLoad = trues(6)
    isFollower = falses(6)
    isTrim = falses(6)

    isLoadOnA = trues(6)
    isFollowerOnA = falses(6)
    isTrimOnA = falses(6)

    isLoadOnb = trues(6)
    isFollowerOnb = falses(6)
    isTrimOnb = falses(6)

    # Initialize flags for BCs having already been specified 
    displacementSpecifiedOnA = falses(6)
    deadLoadSpecifiedOnA = falses(6)
    followerLoadSpecifiedOnA = falses(6)
    displacementSpecifiedOnb = falses(6)
    deadLoadSpecifiedOnb = falses(6)
    followerLoadSpecifiedOnb = falses(6)

    # Initialize generalized displacements and loads as functions of time
    displacementsOnAofTime,deadLoadsOnAofTime,followerLoadsOnAofTime,displacementsOnbofTime,deadLoadsOnbofTime,followerLoadsOnbofTime = Function[t -> 0.0 for _ in 1:6], Function[t -> 0.0 for _ in 1:6], Function[t -> 0.0 for _ in 1:6], Function[t -> 0.0 for _ in 1:6], Function[t -> 0.0 for _ in 1:6], Function[t -> 0.0 for _ in 1:6]

    # Loop generalized displacements/loads indices
    for (i,generalizedDisplacementBasisA,generalizedDisplacementBasisb,deadLoadBasisA,deadLoadBasisb,followerLoadBasisA,followerLoadBasisb) in zip(1:6,[:u1A,:u2A,:u3A,:p1A,:p2A,:p3A],[:u1b,:u2b,:u3b,:p1b,:p2b,:p3b],[:F1A,:F2A,:F3A,:M1A,:M2A,:M3A],[:F1b,:F2b,:F3b,:M1b,:M2b,:M3b],[:Ff1A,:Ff2A,:Ff3A,:Mf1A,:Mf2A,:Mf3A],[:Ff1b,:Ff2b,:Ff3b,:Mf1b,:Mf2b,:Mf3b])

        # Loop specified BCs
        for (type,value,toTrim) in zip(types,values,toBeTrimmed)

            # The specified BC is a generalized displacement in basis A
            if type == generalizedDisplacementBasisA
                # Check if generalized displacement BC was already specified in basis b or if load was already specified in basis A
                if displacementSpecifiedOnb[i]
                    error("Both $generalizedDisplacementBasisA and $generalizedDisplacementBasisb cannot be specified at the same time for $name")
                elseif deadLoadSpecifiedOnA[i]
                    error("Both $generalizedDisplacementBasisA and $deadLoadBasisA cannot be specified at the same time for $name")
                elseif followerLoadSpecifiedOnA[i]
                    error("Both $generalizedDisplacementBasisA and $followerLoadBasisA cannot be specified at the same time for $name")    
                end
                # Update flags
                displacementSpecifiedOnA[i] = true
                isLoadOnA[i] = false
                # Set BC
                displacementsOnAofTime[i] = value
            end

            # The specified BC is a generalized displacement in basis b
            if type == generalizedDisplacementBasisb
                # Check if generalized displacement BC was already specified in basis A
                if displacementSpecifiedOnA[i]
                    error("Both $generalizedDisplacementBasisA and $generalizedDisplacementBasisb cannot be specified at the same time for $name")
                elseif deadLoadSpecifiedOnb[i]
                    error("Both $generalizedDisplacementBasisb and $deadLoadBasisb cannot be specified at the same time for $name")
                elseif followerLoadSpecifiedOnb[i]
                    error("Both $generalizedDisplacementBasisb and $followerLoadBasisb cannot be specified at the same time for $name")     
                end
                # Update flags
                displacementSpecifiedOnb[i] = true
                isLoadOnb[i] = false
                # Set BC
                displacementsOnbofTime[i] = value
            end

            # The specified BC is a dead generalized force on basis A  
            if type == deadLoadBasisA
                # Check if generalized displacement BC was already specified
                if displacementSpecifiedOnA[i]
                    error("Both $generalizedDisplacementBasisA and $deadLoadBasisA cannot be specified at the same time for $name")
                elseif followerLoadSpecifiedOnA[i]
                    error("Please specify $deadLoadBasisA and $followerLoadBasisA in separate BCs")    
                end
                # Update flags
                deadLoadSpecifiedOnA[i] = true 
                # Check if is trim variable
                if toTrim
                    isTrimOnA[i] = true
                end
                # Set BC 
                deadLoadsOnAofTime[i] = value
            end

            # The specified BC is a dead generalized force on basis b
            if type == deadLoadBasisb
                # Check if generalized displacement BC was already specified
                if displacementSpecifiedOnb[i]
                    error("Both $generalizedDisplacementBasisb and $deadLoadBasisb cannot be specified at the same time for $name")
                elseif followerLoadSpecifiedOnb[i]
                    error("Please specify $deadLoadBasisb and $followerLoadBasisb in separate BCs")    
                end
                # Update flags
                deadLoadSpecifiedOnb[i] = true      
                # Check if is trim variable
                if toTrim
                    isTrimOnb[i] = true
                end
                # Set BC
                deadLoadsOnbofTime[i] = value
            end

            # The specified BC is a follower generalized force on basis A  
            if type == followerLoadBasisA
                # Check if generalized displacement BC was already specified
                if displacementSpecifiedOnA[i]
                    error("Both $generalizedDisplacementBasisA and $followerLoadBasisA cannot be specified at the same time for $name")
                elseif deadLoadSpecifiedOnA[i]
                    error("Please specify $deadLoadBasisA and $followerLoadBasisA in separate BCs")    
                end
                # Update flags
                followerLoadSpecifiedOnA[i] = true
                isFollowerOnA[i] = true
                # Check if is trim variable
                if toTrim
                    isTrimOnA[i] = true
                end
                # Set BC
                followerLoadsOnAofTime[i] = value
            end

            # The specified BC is a follower generalized force on basis b  
            if type == followerLoadBasisb
                # Check if generalized displacement BC was already specified
                if displacementSpecifiedOnb[i]
                    error("Both $generalizedDisplacementBasisb and $followerLoadBasisb cannot be specified at the same time for $name")
                elseif deadLoadSpecifiedOnb[i]
                    error("Please specify $deadLoadBasisb and $followerLoadBasisb in separate BCs")    
                end
                # Update TF
                followerLoadSpecifiedOnb[i] = true
                isFollowerOnb[i] = true
                # Check if is trim variable
                if toTrim
                    isTrimOnb[i] = true
                end
                # Set BC
                followerLoadsOnbofTime[i] = value
            end

        end
    end

    # Set flags for being loads, follower loads and trim loads (those already defined on basis A)
    isLoad = isLoadOnA
    isFollower = isFollowerOnA
    isTrim = isTrimOnA

    # Update flags for being loads, follower loads and trim loads on basis A after transfer from basis b
    isLoad = .! (displacementSpecifiedOnA .| vcat(
        [x!=0 for x in R0_n * displacementSpecifiedOnb[1:3]], 
        [x!=0 for x in R0_n * displacementSpecifiedOnb[4:6]]
            ))

    isFollower = isFollower .| vcat(
        [x!=0 for x in R0_n * isFollowerOnb[1:3]], 
        [x!=0 for x in R0_n * isFollowerOnb[4:6]]
            ) 

    isTrim = isTrim .| vcat(
        [x!=0 for x in R0_n * isTrimOnb[1:3]], 
        [x!=0 for x in R0_n * isTrimOnb[4:6]]
            )

    # Initial loads and displacement values
    displacementsOnA = [f(0) for f in displacementsOnAofTime]
    displacementsOnb = [f(0) for f in displacementsOnbofTime]
    deadLoadsOnA = [f(0) for f in deadLoadsOnAofTime]
    deadLoadsOnb = [f(0) for f in deadLoadsOnbofTime]
    followerLoadsOnA = [f(0) for f in followerLoadsOnAofTime]
    followerLoadsOnb = [f(0) for f in followerLoadsOnbofTime]

    # Initialize vector of current BC values
    valuesOnA = displacementsOnA + deadLoadsOnA + followerLoadsOnA
    valuesOnb = displacementsOnb + deadLoadsOnb + followerLoadsOnb
    currentValue = (valuesOnA + vcat(
                    R0_n * valuesOnb[1:3],
                    R0_n * valuesOnb[4:6])
                    ) .* .!isTrim

    # Set deduced data
    deducedData = BCDeducedData(R0_n,isLoad,isFollower,isTrim)

    # Set load functions
    loadFuns = BCLoadFuns(displacementsOnAofTime,deadLoadsOnAofTime,followerLoadsOnAofTime,displacementsOnbofTime,deadLoadsOnbofTime,followerLoadsOnbofTime)
    
    # Initialize state
    state = BCState(displacementsOnA,deadLoadsOnA,followerLoadsOnA,displacementsOnb,deadLoadsOnb,followerLoadsOnb,currentValue)

    return deducedData,loadFuns,state

end


# Updates the boundary conditions' state at the current time index
function update_BC_state!(bc,timeIndex::Int)

    bc.state.displacementsOnA .= bc.history.displacementsOnA[timeIndex]
    bc.state.deadLoadsOnA .= bc.history.deadLoadsOnA[timeIndex]
    bc.state.followerLoadsOnA .= bc.history.followerLoadsOnA[timeIndex]
    bc.state.displacementsOnb .= bc.history.displacementsOnb[timeIndex]
    bc.state.deadLoadsOnb .= bc.history.deadLoadsOnb[timeIndex]
    bc.state.followerLoadsOnb .= bc.history.followerLoadsOnb[timeIndex]
    bc.state.currentValue .= bc.history.currentValue[timeIndex]

end


# Updates the boundary conditions' state at the current time instance
function update_BC_state!(bc,timeNow::Float64)

    # Unpack
    @unpack R0_n,isTrim = bc.deducedData
    @unpack displacementsOnAofTime,deadLoadsOnAofTime,followerLoadsOnAofTime,displacementsOnbofTime,deadLoadsOnbofTime,followerLoadsOnbofTime = bc.loadFuns
    @unpack displacementsOnA,deadLoadsOnA,followerLoadsOnA,displacementsOnb,deadLoadsOnb,followerLoadsOnb,currentValue = bc.state

    # Evaluate prescribed functions
    @inbounds for i in 1:6
        displacementsOnA[i] = displacementsOnAofTime[i](timeNow)
        deadLoadsOnA[i] = deadLoadsOnAofTime[i](timeNow)
        followerLoadsOnA[i] = followerLoadsOnAofTime[i](timeNow)
        displacementsOnb[i] = displacementsOnbofTime[i](timeNow)
        deadLoadsOnb[i] = deadLoadsOnbofTime[i](timeNow)
        followerLoadsOnb[i] = followerLoadsOnbofTime[i](timeNow)
    end

    # Compute current BC values
    valuesOnA = displacementsOnA + deadLoadsOnA + followerLoadsOnA
    valuesOnb = displacementsOnb + deadLoadsOnb + followerLoadsOnb
    currentValue .= (valuesOnA +
                    vcat(
                        R0_n*valuesOnb[1:3],
                        R0_n*valuesOnb[4:6]
                    )) .* .!isTrim

end