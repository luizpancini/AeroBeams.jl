"""
@with_kw mutable struct BC

Boundary conditions composite type

"""
@with_kw mutable struct BC

    # Primary fields (inputs)
    name::String = " "
    beam::Beam 
    node::Int64 
    types::Vector{String} 
    values
    toBeTrimmed::Union{BitVector,Vector{Bool}}

    # Secondary fields (outputs)
    currentValue::Vector{Float64}
    isLoad::BitVector
    isFollower::BitVector
    isTrim::BitVector
    displacementsOnA::Vector{Float64}
    deadLoadsOnA::Vector{Float64}
    followerLoadsOnA::Vector{Float64}
    displacementsOnb::Vector{Float64}
    deadLoadsOnb::Vector{Float64}
    followerLoadsOnb::Vector{Float64}
    initialTrimValue::Vector{Float64}
    R0_n::Matrix{Float64}
    globalNodeID::Int64
    Fmax::Number
    Mmax::Number

end
export BC


"""
create_BC(;name::String="",beam::Beam,node::Int64,types::Vector{String},values,toBeTrimmed::Union{BitVector,Vector{Bool}}=falses(length(types)))

BC constructor

# Keyword arguments
- `name::String` = name of the BC
- `beam::Beam` = beam at which the BC is applied
- `node::Int64` = node of the beam at which the BC is applied
- `types::Vector{String} `= types of BCs applied to the node (generalized forces and displacements)
- `values` = corresponding values of the applied BCs (constants or functions of time)
- `toBeTrimmed::Union{BitVector,Vector{Bool}}` = TF on whether the BC is to be trimmed
"""
function create_BC(;name::String="",beam::Beam,node::Int64,types::Vector{String},values,toBeTrimmed::Union{BitVector,Vector{Bool}}=falses(length(types)))

    # Validate inputs
    @assert 1 <= node <= beam.nElements+1 "the BC'ed beam does not contain the BC'ed node"
    @assert types == unique!(deepcopy(types)) "there are repeated BC types"
    @assert length(types) == length(values) "set one value for each BC type"
    @assert length(types) == length(toBeTrimmed) "set a toBeTrimmed value for each BC type"
    for (i,type) in enumerate(types)
        @assert in(type,["u1A","u2A","u3A","p1A","p2A","p3A","u1b","u2b","u3b","p1b","p2b","p3b","F1A","F2A","F3A","M1A","M2A","M3A","F1b","F2b","F3b","M1b","M2b","M3b","Ff1A","Ff2A","Ff3A","Mf1A","Mf2A","Mf3A","Ff1b","Ff2b","Ff3b","Mf1b","Mf2b","Mf3b"])
        if !(type in ["F1A","F2A","F3A","M1A","M2A","M3A","F1b","F2b","F3b","M1b","M2b","M3b","Ff1A","Ff2A","Ff3A","Mf1A","Mf2A","Mf3A","Ff1b","Ff2b","Ff3b","Mf1b","Mf2b","Mf3b"]) && toBeTrimmed[i]
            error("Trim variables must be loads")
        end
    end
    for value in values
        if !(value isa Number || value isa Function)
            throw(ArgumentError("Elements of 'values' must be either numbers or functions"))
        end
    end
    toBeTrimmed = BitVector(toBeTrimmed)

    # Get TFs for types of BC and validate
    containsDisplacementA = any([type in ["u1A","u2A","u3A","p1A","p2A","p3A"] for type in types])
    containsDisplacementb = any([type in ["u1b","u2b","u3b","p1b","p2b","p3b"] for type in types])
    containsDeadLoadA = any([type in ["F1A","F2A","F3A","M1A","M2A","M3A"] for type in types])
    containsDeadLoadb = any([type in ["F1b","F2b","F3b","M1b","M2b","M3b"] for type in types])
    containsFollowerLoadA = any([type in ["Ff1A","Ff2A","Ff3A","Mf1A","Mf2A","M3A"] for type in types])
    containsFollowerLoadb = any([type in ["Ff1b","Ff2b","Ff3b","Mf1b","Mf2b","M3b"] for type in types])
    if containsDisplacementA && containsDisplacementb
        error("Specify displacements and rotations either solely on basis A or solely on basis b")
    end
    if containsDeadLoadA && containsDeadLoadb
        error("Please input dead loads on bases A and b on separate BCs")
    end
    if containsFollowerLoadA && containsFollowerLoadb
        error("Please input follower loads on bases A and b on separate BCs")
    end

    # Get nodal rotation tensor from basis A to basis b
    R0_n = node == length(beam.elements)+1 ? beam.elements[node-1].R0_n2 : beam.elements[node].R0_n1

    # globalNodeID is updated upon assembly of the model

    # Initialize BC
    bc = BC(name,beam,node,types,values,toBeTrimmed,Vector{Float64}(),BitVector(),BitVector(),BitVector(),Vector{Float64}(),Vector{Float64}(),Vector{Float64}(),Vector{Float64}(),Vector{Float64}(),Vector{Float64}(),Vector{Float64}(),R0_n,0,0,0)

    # Get initial BC data
    update_BC_data!(bc)

    return bc
end
export create_BC


"""
update_BC_data!(bc::BC,timeNow::Number=0)

Updates the boundary conditions at the current time

"""
function update_BC_data!(bc::BC,timeNow::Number=0)

    @unpack name,types,values,R0_n,toBeTrimmed,Fmax,Mmax = bc

    # Initialize outputs (default corresponds to no BCs applied on basis A)
    currentValue = zeros(6)
    isLoad = trues(6)
    isFollower = falses(6)
    isTrim = falses(6)

    # Initialize vectors in basis A
    displacementsOnA = zeros(6)
    deadLoadsOnA = zeros(6)
    followerLoadsOnA = zeros(6)
    isLoadOnA = trues(6)
    isFollowerOnA = falses(6)
    isTrimOnA = falses(6)

    # Initialize vectors in basis b
    displacementsOnb = zeros(6)
    deadLoadsOnb = zeros(6)
    followerLoadsOnb = zeros(6)
    isLoadOnb = trues(6)
    isFollowerOnb = falses(6)
    isTrimOnb = falses(6)

    # Initialize TF for BCs having already being specified 
    displacementSpecifiedOnA = falses(6)
    deadLoadSpecifiedOnA = falses(6)
    followerLoadSpecifiedOnA = falses(6)
    displacementSpecifiedOnb = falses(6)
    deadLoadSpecifiedOnb = falses(6)
    followerLoadSpecifiedOnb = falses(6)

    # Loop generalized displacements/loads indices
    for (i,generalizedDisplacementBasisA,generalizedDisplacementBasisb,deadLoadBasisA,deadLoadBasisb,followerLoadBasisA,followerLoadBasisb) in zip(1:6,["u1A","u2A","u3A","p1A","p2A","p3A"],["u1b","u2b","u3b","p1b","p2b","p3b"],["F1A","F2A","F3A","M1A","M2A","M3A"],["F1b","F2b","F3b","M1b","M2b","M3b"],["Ff1A","Ff2A","Ff3A","Mf1A","Mf2A","Mf3A"],["Ff1b","Ff2b","Ff3b","Mf1b","Mf2b","Mf3b"])

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
                # Update TFs
                displacementSpecifiedOnA[i] = true
                isLoadOnA[i] = false
                # Check whether the BC is a constant or a function of time and set current value accordingly
                displacementsOnA[i] = value isa Function ? value(timeNow) : value     
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
                # Update TFs
                displacementSpecifiedOnb[i] = true
                isLoadOnb[i] = false
                # Check whether the BC is a constant or a function of time and set current value accordingly
                displacementsOnb[i] = value isa Function ? value(timeNow) : value      
            end

            # The specified BC is a dead generalized force on basis A  
            if type == deadLoadBasisA
                # Check if generalized displacement BC was already specified
                if displacementSpecifiedOnA[i]
                    error("Both $generalizedDisplacementBasisA and $deadLoadBasisA cannot be specified at the same time for $name")
                elseif followerLoadSpecifiedOnA[i]
                    error("Please specify $deadLoadBasisA and $followerLoadBasisA in separate BCs")    
                end
                # Update TFs
                deadLoadSpecifiedOnA[i] = true
                # Check whether the BC is a constant or a function of time and set current value accordingly
                deadLoadsOnA[i] = value isa Function ? value(timeNow) : value   
                # Check if is trim variable
                if toTrim
                    isTrimOnA[i] = true
                end
            end

            # The specified BC is a dead generalized force on basis b
            if type == deadLoadBasisb
                # Check if generalized displacement BC was already specified
                if displacementSpecifiedOnb[i]
                    error("Both $generalizedDisplacementBasisb and $deadLoadBasisb cannot be specified at the same time for $name")
                elseif followerLoadSpecifiedOnb[i]
                    error("Please specify $deadLoadBasisb and $followerLoadBasisb in separate BCs")    
                end
                # Update TFs
                deadLoadSpecifiedOnb[i] = true
                # Check whether the BC is a constant or a function of time and set current value accordingly
                deadLoadsOnb[i] = value isa Function ? value(timeNow) : value         
                # Check if is trim variable
                if toTrim
                    isTrimOnb[i] = true
                end
            end

            # The specified BC is a follower generalized force on basis A  
            if type == followerLoadBasisA
                # Check if generalized displacement BC was already specified
                if displacementSpecifiedOnA[i]
                    error("Both $generalizedDisplacementBasisA and $followerLoadBasisA cannot be specified at the same time for $name")
                elseif deadLoadSpecifiedOnA[i]
                    error("Please specify $deadLoadBasisA and $followerLoadBasisA in separate BCs")    
                end
                # Update TFs
                followerLoadSpecifiedOnA[i] = true
                isFollowerOnA[i] = true
                # Check whether the BC is a constant or a function of time and set current value accordingly
                followerLoadsOnA[i] = value isa Function ? value(timeNow) : value
                # Check if is trim variable
                if toTrim
                    isTrimOnA[i] = true
                end
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
                # Check whether the BC is a constant or a function of time and set current value accordingly
                followerLoadsOnb[i] = value isa Function ? value(timeNow) : value     
                # Check if is trim variable
                if toTrim
                    isTrimOnb[i] = true
                end
            end

        end
    end

    # Set current TFs for being loads (those already defined on basis A)
    isLoad = isLoadOnA

    # Set current TFs for being follower loads (those already defined on basis A)
    isFollower = isFollowerOnA

    # Set current TFs for being trim loads (those already defined on basis A)
    isTrim = isTrimOnA

    # Set current prescribed values (as those already defined on basis A, excluding trim values)
    currentValue = (displacementsOnA + deadLoadsOnA + followerLoadsOnA) .* .!isTrim

    # Update current TFs for being loads on basis A (after transfer from basis b)
    isLoad = isLoad .& vcat([x!=0 for x in R0_n * isLoadOnb[1:3]], [x!=0 for x in R0_n * isLoadOnb[4:6]]) 

    # Update current TFs for being follower loads on basis A (after transfer from basis b)
    isFollower = isFollower .| vcat([x!=0 for x in R0_n * isFollowerOnb[1:3]], [x!=0 for x in R0_n * isFollowerOnb[4:6]]) 

    # Update current TFs for being trim loads on basis A (after transfer from basis b)
    isTrim = isTrim .| vcat([x!=0 for x in R0_n * isTrimOnb[1:3]], [x!=0 for x in R0_n * isTrimOnb[4:6]])

    # Update current prescribed values (excludes trim values) on basis A (after transfer from basis b)
    currentValue += vcat(R0_n * (displacementsOnb[1:3]+deadLoadsOnb[1:3]+followerLoadsOnb[1:3]), R0_n * (displacementsOnb[4:6]+deadLoadsOnb[4:6]+followerLoadsOnb[4:6])) .* .!isTrim

    # Update maximum values of forces and moments
    Fmax = max(Fmax,maximum(abs.(currentValue[1:3])))
    Mmax = max(Mmax,maximum(abs.(currentValue[4:6])))

    # Pack
    @pack! bc = currentValue,isLoad,isFollower,isTrim,displacementsOnA,deadLoadsOnA,followerLoadsOnA,displacementsOnb,deadLoadsOnb,followerLoadsOnb,Fmax,Mmax

    # If in initial time
    if timeNow == 0
        # Set initial trim values (those already defined on basis A)
        initialTrimValue = (deadLoadsOnA + followerLoadsOnA) .* isTrim
        
        # Update initial trim values on basis A
        initialTrimValue += vcat(R0_n * (deadLoadsOnb[1:3]+followerLoadsOnb[1:3]), R0_n * (deadLoadsOnb[4:6]+followerLoadsOnb[4:6])) .* isTrim
        
        @pack! bc = initialTrimValue
    end

end

