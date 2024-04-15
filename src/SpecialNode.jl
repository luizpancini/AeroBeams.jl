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
    connectedElementsID::Vector{Int64} 
    connectedElements::Vector{Element} 
    ζonElements::Vector{Int64} 
    BCs::Vector{BC} = Vector{BC}()
    isLoad::Vector{Bool} = trues(6)
    eqs_Fu::Vector{Vector{Int64}} = Vector{Vector{Int64}}()
    eqs_Fp::Vector{Vector{Int64}} = Vector{Vector{Int64}}()
    eqs_FF::Vector{Vector{Int64}} = Vector{Vector{Int64}}()
    eqs_FM::Vector{Vector{Int64}} = Vector{Vector{Int64}}()
    eqs_FF_sep::Vector{Vector{Int64}} = Vector{Vector{Int64}}() 
    eqs_FM_sep::Vector{Vector{Int64}} = Vector{Vector{Int64}}()
    DOF_uF::Vector{Int64} = Vector{Int64}()
    DOF_pM::Vector{Int64} = Vector{Int64}()
    DOF_trimLoads::Vector{Int64} = zeros(Int64,6)

end

# Constructor for special nodes with BCs
function SpecialNode(localID::Int64,globalID::Int64,connectedElementsID::Vector{Int64},connectedElements::Vector{Element},ζonElements::Vector{Int64},BCs::Vector{BC})

    # Set isLoad TF
    isLoad = trues(6)
    for BC in BCs
        isLoad = isLoad .& BC.isLoad
    end

    return SpecialNode(localID,globalID,connectedElementsID,connectedElements,ζonElements,BCs,isLoad,[Vector{Int64}() for _ in 1:length(connectedElements)],[Vector{Int64}() for _ in 1:length(connectedElements)],[Vector{Int64}() for _ in 1:length(connectedElements)],[Vector{Int64}() for _ in 1:length(connectedElements)],[Vector{Int64}() for _ in 1:length(connectedElements)],[Vector{Int64}() for _ in 1:length(connectedElements)],Vector{Int64}(),Vector{Int64}(),zeros(Int64,6))

end

# Constructor for special nodes without BCs
function SpecialNode(localID::Int64,globalID::Int64,connectedElementsID::Vector{Int64},connectedElements::Vector{Element},ζonElements::Vector{Int64})

    return SpecialNode(localID,globalID,connectedElementsID,connectedElements,ζonElements,Vector{BC}(),trues(6),[Vector{Int64}() for _ in 1:length(connectedElements)],[Vector{Int64}() for _ in 1:length(connectedElements)],[Vector{Int64}() for _ in 1:length(connectedElements)],[Vector{Int64}() for _ in 1:length(connectedElements)],[Vector{Int64}() for _ in 1:length(connectedElements)],[Vector{Int64}() for _ in 1:length(connectedElements)],Vector{Int64}(),Vector{Int64}(),zeros(Int64,6))

end