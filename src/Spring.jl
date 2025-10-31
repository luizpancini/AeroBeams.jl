#
# @with_kw mutable struct Spring

#     Spring composite type

#
@with_kw mutable struct Spring

    # Primary (inputs to spring creation)
    basis::String
    elementsIDs::Vector{Int}
    nodesSides::Vector{Int}
    ku::Vector{<:Real}
    kp::Vector{<:Real}

    # Secondary (outputs from spring creation)
    hasDoubleAttachment::Bool
    nodesSpecialIDs::Vector{Int}
    nodesGlobalIDs::Vector{Int}
    Ku::Matrix{Float64}
    Kp::Matrix{Float64}
    Fs::Vector{Float64}
    Ms::Vector{Float64}
    Δu::Vector{Float64}
    Δp::Vector{Float64}

end
export Spring


"""
    create_Spring(; kwargs...)

Creates a spring

# Keyword arguments
- `basis::String`: basis on which stiffnesses are defined ("b" or "A")
- `elementsIDs::Vector{Int}`: local IDs of the element(s)' node(s) to which the spring is attached
- `nodesSides::Vector{Int}`: sides (1 or 2) of the node(s) to which the spring is attached
- `ku::Vector{<:Real}`: translational stiffness vector, resolved in basis A
- `kp::Vector{<:Real}`: rotational stiffness vector, resolved in basis A
"""
function create_Spring(; basis::String="A",elementsIDs::Vector{Int},nodesSides::Vector{Int},ku::Vector{<:Real}=zeros(3),kp::Vector{<:Real}=zeros(3))

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
    @assert all(x -> x>=0, [ku; kp]) "elements of ku and kp must be positive"

    # TF for being attached to beams at both ends
    hasDoubleAttachment = N == 2

    # Initialize IDs of attachment nodes on list of special nodes (updated later on model assembly)
    nodesSpecialIDs = Vector{Int}()

    # Initialize global IDs of attachment nodes (updated later on model assembly)
    nodesGlobalIDs = Vector{Int}()

    # Initialize spring stiffness matrices, resolved in basis A
    Ku = [ku[1] 0 0; 0 ku[2] 0; 0 0 ku[3]]
    Kp = [kp[1] 0 0; 0 kp[2] 0; 0 0 kp[3]]

    # Initialize generalized spring displacements, resolved in basis A
    Δu = Δp = zeros(3)

    # Initialize spring load vectors, resolved in basis A
    Fs = Ms = zeros(3)

    return Spring(basis=basis,elementsIDs=elementsIDs,nodesSides=nodesSides,ku=ku,kp=kp,hasDoubleAttachment=hasDoubleAttachment,nodesSpecialIDs=nodesSpecialIDs,nodesGlobalIDs=nodesGlobalIDs,Ku=Ku,Kp=Kp,Fs=Fs,Ms=Ms,Δu=Δu,Δp=Δp)

end
export create_Spring