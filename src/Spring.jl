"""
@with_kw mutable struct Spring

    Spring composite type

# Fields
- 
"""
@with_kw mutable struct Spring

    # Primary (inputs to spring creation)
    basis::String
    elementsIDs::Vector{Int64}
    nodesSides::Vector{Int64}
    ku::Vector{<:Number}
    kp::Vector{<:Number}
    kTranslational::Number
    kTwist::Number
    kIPBending::Number
    kOOPBending::Number

    # Secondary (outputs from spring creation)
    hasDoubleAttachment::Bool
    nodesSpecialIDs::Vector{Int64}
    nodesGlobalIDs::Vector{Int64}
    Ku::Matrix{Float64}
    Kp::Matrix{Float64}

end

# Constructor
function create_Spring(; basis::String="A",elementsIDs::Vector{Int64},nodesSides::Vector{Int64},ku::Vector{<:Number}=zeros(3),kp::Vector{<:Number}=zeros(3),kTranslational::Number=0,kTwist::Number=0,kIPBending::Number=0,kOOPBending::Number=0)

    # Number of attachments
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
        Kp = [kTwist 0 0; 0 kIPBending 0; 0 0 kOOPBending]
    else
        Ku = [ku[1] 0 0; 0 ku[2] 0; 0 0 ku[3]]
        Kp = [kp[1] 0 0; 0 kp[2] 0; 0 0 kp[3]]
    end

    return Spring(basis=basis,elementsIDs=elementsIDs,nodesSides=nodesSides,ku=ku,kp=kp,kTranslational=kTranslational,kTwist=kTwist,kIPBending=kIPBending,kOOPBending=kOOPBending,hasDoubleAttachment=hasDoubleAttachment,nodesSpecialIDs=nodesSpecialIDs,nodesGlobalIDs=nodesGlobalIDs,Ku=Ku,Kp=Kp)

end
export create_Spring