#
# @with_kw mutable struct Spring

#     Spring composite type

#
@with_kw mutable struct Spring

    # Primary (inputs to spring creation)
    basis::String
    elementsIDs::Vector{Int64}
    nodesSides::Vector{Int64}
    ku::Vector{<:Real}
    kp::Vector{<:Real}
    kTranslational::Real
    kTwist::Real
    kIPBending::Real
    kOOPBending::Real

    # Secondary (outputs from spring creation)
    hasDoubleAttachment::Bool
    nodesSpecialIDs::Vector{Int64}
    nodesGlobalIDs::Vector{Int64}
    Ku::Matrix{Float64}
    Kp::Matrix{Float64}

end


"""
    create_Spring(; kwargs...)

Creates a spring

# Keyword arguments
- `basis::String` = basis on which stiffnesses are defined
- `elementsIDs::Vector{Int64}` = local IDs of the element(s)' node(s) to which the spring is attached
- `nodesSides::Vector{Int64}` = sides (1 or 2) of the node(s) to which the spring is attached
- `ku::Vector{<:Real}` = translational stiffness vector
- `kp::Vector{<:Real}` = rotational stiffness vector
- `kTranslational::Real` = translational stiffness (same in all directions)
- `kTwist::Real` = twist stiffness (rotation about x1 direction)
- `kIPBending::Real` = in-plane bending stiffness (rotation about x3 direction)
- `kOOPBending::Real` = out-of-plane bending stiffness (rotation about x2 direction)
"""
function create_Spring(; basis::String="A",elementsIDs::Vector{Int64},nodesSides::Vector{Int64},ku::Vector{<:Real}=zeros(3),kp::Vector{<:Real}=zeros(3),kTranslational::Real=0,kTwist::Real=0,kIPBending::Real=0,kOOPBending::Real=0)

    # Real of attachments
    N = length(elementsIDs)

    # Validate
    @assert basis in ["b","A"] "basis must be either 'b' or 'A'"
    @assert length(elementsIDs) == length(nodesSides) "assign a node side for each element"
    @assert 1 <= N <= 2 "spring must be attached to at least one and at most two nodes"
    for (elementID,nodeSide) in zip(elementsIDs,nodesSides) 
        @assert elementID > 0 "elementID must be positive"
        @assert nodeSide in [1,2] "nodeSide must be either 1 or 2"
    end
    @assert length(ku) == length(kp) == 3 "ku and kp must be three-element vectors"
    @assert all(x -> x>=0, [ku; kp]) "elements of ku and kp must be >= 0"
    @assert all(x -> x>=0, [kTranslational; kTwist; kIPBending; kOOPBending]) "kTranslational, kTwist, kIPBending and kOOPBending must be >= 0"

    # TF for being attached to beams at both ends
    hasDoubleAttachment = N == 2

    # Initialize IDs of attachment nodes on list of special nodes (updated later on model assembly)
    nodesSpecialIDs = Vector{Int64}()

    # Initialize global IDs of attachment nodes (updated later on model assembly)
    nodesGlobalIDs = Vector{Int64}()

    # Initialize spring stiffness matrices, resolved in basis A
    if hasDoubleAttachment
        Ku = [kTranslational 0 0; 0 kTranslational 0; 0 0 kTranslational]
        Kp = [kTwist 0 0; 0 kOOPBending 0; 0 0 kIPBending]
    else
        Ku = [ku[1] 0 0; 0 ku[2] 0; 0 0 ku[3]]
        Kp = [kp[1] 0 0; 0 kp[2] 0; 0 0 kp[3]]
    end

    return Spring(basis=basis,elementsIDs=elementsIDs,nodesSides=nodesSides,ku=ku,kp=kp,kTranslational=kTranslational,kTwist=kTwist,kIPBending=kIPBending,kOOPBending=kOOPBending,hasDoubleAttachment=hasDoubleAttachment,nodesSpecialIDs=nodesSpecialIDs,nodesGlobalIDs=nodesGlobalIDs,Ku=Ku,Kp=Kp)

end
export create_Spring