"""
@with_kw mutable struct SpecialNode

    Special node composite type

# Fields
- 
"""
@with_kw mutable struct SpecialNode

    # Fields
    localID::Int64
    globalID::Int64 
    connectedElementsGlobalIDs::Vector{Int64} 
    connectedElements::Vector{Element} 
    ζonElements::Vector{Int64} 
    BCs::Vector{BC}
    uIsPrescribed::BitVector
    pIsPrescribed::BitVector
    trimIsPrescribed::BitVector
    eqs_Fu::Vector{Vector{Int64}} = [Vector{Int64}() for _ in 1:length(connectedElements)]
    eqs_Fp::Vector{Vector{Int64}} = [Vector{Int64}() for _ in 1:length(connectedElements)]
    eqs_FF::Vector{Vector{Int64}} = [Vector{Int64}() for _ in 1:length(connectedElements)]
    eqs_FM::Vector{Vector{Int64}} = [Vector{Int64}() for _ in 1:length(connectedElements)]
    eqs_FF_sep::Vector{Vector{Int64}} = [Vector{Int64}() for _ in 1:length(connectedElements)]
    eqs_FM_sep::Vector{Vector{Int64}} = [Vector{Int64}() for _ in 1:length(connectedElements)]
    DOF_uF::Vector{Int64} = Vector{Int64}()
    DOF_pM::Vector{Int64} = Vector{Int64}()
    DOF_trimLoads::Vector{Int64} = Vector{Int64}()
    u::Vector{Float64} = Vector{Float64}()
    p::Vector{Float64} = Vector{Float64}()
    F::Vector{Float64} = Vector{Float64}()
    M::Vector{Float64} = Vector{Float64}()
    F_p::Matrix{Float64} = zeros(3,3)
    M_p::Matrix{Float64} = zeros(3,3)

end

# Constructor 
function SpecialNode(localID::Int64,globalID::Int64,connectedElementsGlobalIDs::Vector{Int64},connectedElements::Vector{Element},ζonElements::Vector{Int64},BCs::Vector{BC}=Vector{BC}())

    # Set TFs for generalized displacements and trim loads being prescribed 
    isLoad = trues(6)
    isTrim = falses(6)
    for BC in BCs
        isLoad = isLoad .& BC.isLoad
        isTrim = isTrim .| BC.isTrim
    end
    uIsPrescribed = .!isLoad[1:3]
    pIsPrescribed = .!isLoad[4:6]
    trimIsPrescribed = isTrim

    return SpecialNode(localID=localID,globalID=globalID,connectedElementsGlobalIDs=connectedElementsGlobalIDs,connectedElements=connectedElements,ζonElements=ζonElements,BCs=BCs,uIsPrescribed=uIsPrescribed,pIsPrescribed=pIsPrescribed,trimIsPrescribed=trimIsPrescribed)

end