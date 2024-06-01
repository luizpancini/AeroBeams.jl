"""
@with_kw mutable struct Spring

    Spring composite type

# Fields
- 
"""
@with_kw mutable struct Spring

    # Primary (inputs to spring creation)
    elementID::Int64
    localNode::Int64
    ku::Vector{<:Number}
    kp::Vector{<:Number}

    # Secondary (outputs from spring creation)
    nodeGlobalID::Int64
    R0::Matrix{Float64}

end

# Constructor
function create_Spring(; elementID::Int64,localNode::Int64,ku::Vector{<:Number}=zeros(3),kp::Vector{<:Number}=zeros(3))

    # Validate
    @assert elementID > 0 "elementID must be positive"
    @assert localNode in [1,2] "localNode must be either 1 or 2"
    @assert all(x -> x>=0, ku) "elements of ku must be >= 0"
    @assert all(x -> x>=0, kp) "elements of kp must be >= 0"

    # Initialize node global ID (updated later on model assembly)
    nodeGlobalID = 0

    # Initialize nodal rotation tensor from basis A to basis b
    R0 = I3

    return Spring(elementID=elementID,localNode=localNode,ku=ku,kp=kp,nodeGlobalID=nodeGlobalID,R0=R0)

end
export create_Spring