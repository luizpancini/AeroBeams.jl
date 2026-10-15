## SPRING COMPOSITE TYPES
# ------------------------------------------------------------------------------
# Function types
const SpringForceDefFun = FunctionWrappers.FunctionWrapper{SVector{3,Float64}, Tuple{SMatrix{3,3,Float64,9},SVector{3,Float64}}}
const SpringForceFun = FunctionWrappers.FunctionWrapper{SVector{3,Float64}, Tuple{SVector{3,Float64}}}
const SpringJacFun = FunctionWrappers.FunctionWrapper{SMatrix{3,3,Float64,9}, Tuple{SVector{3,Float64}}}

# Spring constitutive data composite type
# Fields
# - `ku`: translational stiffness constant
# - `kp`: rotational stiffness constant
# - `forceFunction`: function of two inputs (stiffness constant and deformation) that outputs the spring force
# - `momentFunction`: function of two inputs (stiffness constant and deformation) that outputs the spring moment
#
struct SpringConstitutiveData

    ku::Float64
    kp::Float64

    forceFunction::SpringForceDefFun
    momentFunction::SpringForceDefFun

end


# Spring assembly data composite type
# Fields
# - `Δu0`: translational pre-deformation
# - `Δp0`: rotational pre-deformation
# - `hasDoubleAttachment`: flag for being attached to beams at both ends
# - `beamConstructionIDs`: construction IDs of the beams to which the spring is attached
# - `nodesLocalIDs`: local IDs (on the beams) of the node(s) to which the spring is attached
# - `elemsLocalIDs`: beam element(s) local ID(s) to which the spring is attached
# - `nodesSides`: side of the spring's attachment node(s) on the beam element(s)
# - `nodesGlobalIDs`: global IDs of the spring's attached nodes
# - `nodesSpecialIDs`: special node IDs of the spring's attached nodes
#
struct SpringAssemblyData

    Δu0::Float64
    Δp0::Float64

    hasDoubleAttachment::Bool

    beamConstructionIDs::Vector{Int}

    nodesLocalIDs::Vector{Int}
    elemsLocalIDs::Vector{Int}
    nodesSides::Vector{Int}

    nodesGlobalIDs::Vector{Int}
    nodesSpecialIDs::Vector{Int}

end
SpringAssemblyData() = SpringAssemblyData(
    0.0,
    0.0,
    false,
    Vector{Int}(),
    ntuple(_ -> Int[], 3)...,
    ntuple(_ -> Int[], 2)...,
    )

SpringAssemblyData(spring, hasDoubleAttachment::Bool, beamConstructionIDs::Vector{Int}, nodesLocalIDs::AbstractVector{Int}, elemsLocalIDs::AbstractVector{Int}, nodesSides::AbstractVector{Int}) =
    SpringAssemblyData(
        spring.assemblyData.Δu0, spring.assemblyData.Δp0,
        hasDoubleAttachment,
        beamConstructionIDs,
        collect(nodesLocalIDs),collect(elemsLocalIDs),collect(nodesSides),
        spring.assemblyData.nodesGlobalIDs,spring.assemblyData.nodesSpecialIDs,
    )

SpringAssemblyData(spring, nodesGlobalIDs::AbstractVector{Int}) =
    SpringAssemblyData(
        spring.assemblyData.Δu0, spring.assemblyData.Δp0,
        spring.assemblyData.hasDoubleAttachment,
        spring.assemblyData.beamConstructionIDs,
        spring.assemblyData.nodesLocalIDs, spring.assemblyData.elemsLocalIDs, spring.assemblyData.nodesSides,
        collect(nodesGlobalIDs),spring.assemblyData.nodesSpecialIDs,
    )

SpringAssemblyData(spring, nodesGlobalIDs::AbstractVector{Int},nodesSpecialIDs::AbstractVector{Int}) =
    SpringAssemblyData(
        spring.assemblyData.Δu0, spring.assemblyData.Δp0,
        spring.assemblyData.hasDoubleAttachment,
        spring.assemblyData.beamConstructionIDs,
        spring.assemblyData.nodesLocalIDs, spring.assemblyData.elemsLocalIDs, spring.assemblyData.nodesSides,
        collect(nodesGlobalIDs),collect(nodesSpecialIDs),
    )


# Spring assembled constitutive data composite type
# Fields
# - `Ku`: translational stiffness matrix, resolved in basis A
# - `Kp`: rotational stiffness matrix, resolved in basis A
# - `FsFun`: function of spring translational deformation vector that outputs the spring force vector resolved in basis A
# - `MsFun`: function of spring rotational deformation vector that outputs the spring moment vector resolved in basis A
# - `Fs_ΔuFun`: function of spring translational deformation vector that outputs the derivative of the spring force vector resolved in basis A
# - `Ms_ΔpFun`: function of spring rotational deformation vector that outputs the derivative of the spring moment vector resolved in basis A
#
struct SpringAssembledConstitutiveData

    Ku::SMatrix{3,3,Float64,9}
    Kp::SMatrix{3,3,Float64,9}

    FsFun::SpringForceFun
    MsFun::SpringForceFun
    Fs_ΔuFun::SpringJacFun
    Ms_ΔpFun::SpringJacFun

end
SpringAssembledConstitutiveData() = begin
    zeroFun = x -> zeros(SVector{3,Float64})
    SpringAssembledConstitutiveData(
        zeros(SMatrix{3,3,Float64,9}), zeros(SMatrix{3,3,Float64,9}),
        ntuple(_ -> zeroFun, 4)...,
    )
end


# Spring state composite type
# Fields
# - `Fs`: spring force vector, resolved in basis A
# - `Ms`: spring moment vector, resolved in basis A
# - `Δu`: spring translational deformation vector, resolved in basis A
# - `Δp`: spring rotational deformation vector, resolved in basis A
#
struct SpringState

    Fs::MVector{3,Float64}
    Ms::MVector{3,Float64}
    Δu::MVector{3,Float64}
    Δp::MVector{3,Float64}

end
SpringState() = SpringState(
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
)


# DoubleAttachmentSpringData composite type
# Fields
# - `otherNodeSpecialID`: special ID of the other node attached to the spring
# - `otherNodeDOF_uF`: global DOFs of the translational displacements and forces of the other node
# - `otherNodeDOF_pM`: global DOFs of the rotational displacements and moments of the other node
# - `uOtherNode`: translational displacement vector of the other node
# - `pOtherNode`: rotation parameters vector of the other node
#
struct DoubleAttachmentSpringData

    otherNodeSpecialID::Int
    otherNodeDOF_uF::SVector{3,Int}
    otherNodeDOF_pM::SVector{3,Int}

    uOtherNode::MVector{3,Float64}
    pOtherNode::MVector{3,Float64}

end
DoubleAttachmentSpringData(otherNodeSpecialID::Int, otherNodeDOF_uF::SVector{3,Int}, otherNodeDOF_pM::SVector{3,Int}) =
    DoubleAttachmentSpringData(
        otherNodeSpecialID, otherNodeDOF_uF, otherNodeDOF_pM,
        ntuple(_ -> zeros(MVector{3,Float64}), 2)..., 
    )


# Spring composite type
# Fields
# - `constitutive`: constitutive data of the spring
# - `assemblyData`: assembly data for the spring
# - `assemblyConstitutive`: assembled constitutive data for the spring
# - `state`: state of the spring
#
@with_kw mutable struct Spring

    constitutive::SpringConstitutiveData

    assemblyData::SpringAssemblyData = SpringAssemblyData()

    assemblyConstitutive::SpringAssembledConstitutiveData = SpringAssembledConstitutiveData()

    state::SpringState = SpringState()

end
export Spring


## PUBLIC API
# ------------------------------------------------------------------------------
"""
    create_Spring(; kwargs...)

Creates a spring, whose displacement is assumed to be small (changes in direction of the spring's loads are neglected)

# Keyword arguments
- `ku::Real`: translational stiffness constant
- `kp::Real`: rotational stiffness constant
- `forceFunction::Union{Missing,<:Function}`: function of two inputs (stiffness constant and deformation) that outputs the spring force. Default is a linear function: `f(k,x) = k * x` where `k=ku`
- `momentFunction::Union{Missing,<:Function}`: function of two inputs (stiffness constant and deformation) that outputs the spring moment. Default is a linear function: `f(k,x) = k * x` where `k=kp`
"""
function create_Spring(;
    ku::Real=0,
    kp::Real=0,
    forceFunction::Union{Missing,<:Function}=missing,
    momentFunction::Union{Missing,<:Function}=missing,
    )

    # Validate and set correct types
    @assert ku >= 0 "ku must be non-negative"
    @assert kp >= 0 "kp must be non-negative"
    if !ismissing(forceFunction)
        @assert applicable(forceFunction, 0, 0) "forceFunction must be callable with two arguments (k,x)"
        @assert length(forceFunction(randn(3,3),randn(3))) == 3 "forceFunction(k,x) must return a 3-element vector given k as a 3x3 stiffness matrix and x as a 3-element deformation vector"
    end
    if !ismissing(momentFunction)
        @assert applicable(momentFunction, 0, 0) "momentFunction must be callable with two arguments (k,x)"
        @assert length(momentFunction(randn(3,3),randn(3))) == 3 "momentFunction(k,x) must return a 3-element vector given k as a 3x3 stiffness matrix and x as a 3-element deformation vector"
    end
    ku = Float64(ku)
    kp = Float64(kp)

    # Set default (linear) force and moment functions if not provided
    if ismissing(forceFunction)
        forceFunction = (k, x) -> k * x
    end
    if ismissing(momentFunction)
        momentFunction = (k, x) -> k * x
    end

    # Wrap into FunctionWrappers to concretize the field types
    forceFunctionWrapped = wrap_spring_def_fun(forceFunction)
    momentFunctionWrapped = wrap_spring_def_fun(momentFunction)

    # Initialize spring constitutive data
    constitutive = SpringConstitutiveData(ku,kp,forceFunctionWrapped,momentFunctionWrapped)

    return Spring(
        constitutive=constitutive,
        )

end
export create_Spring


## PRIVATE API
# ------------------------------------------------------------------------------
# Function barrier: concretizes f's type before it's captured by the FunctionWrapper's closure (avoids a residual dynamic dispatch inside the wrapper)
@inline wrap_spring_def_fun(f::F) where F = 
    SpringForceDefFun((K,x) -> SVector{3,Float64}(f(K,x)))

# Non-allocating central finite-difference Jacobian for spring loads w.r.t. deformations
@inline function jacobian3x3_fd(f::F, x::SVector{3,Float64}) where F
    h = 1e-6
    c1 = (f(x + SVector(h,0,0)) - f(x - SVector(h,0,0))) / 2h
    c2 = (f(x + SVector(0,h,0)) - f(x - SVector(0,h,0))) / 2h
    c3 = (f(x + SVector(0,0,h)) - f(x - SVector(0,0,h))) / 2h
    return hcat(c1,c2,c3)
end

# Updates assembly and constitutive data of a spring
function update_spring_assembly_data!(beams,spring,nodesLocalIDs,Δu0,Δp0,p=z3,basis=:b)

    @unpack ku,kp,forceFunction,momentFunction = spring.constitutive

    # Number of attachments
    N = length(nodesLocalIDs)

    # Flag for being attached to beams at both ends
    hasDoubleAttachment = N == 2

    # Array of local IDs of the first beam elements containing the spring's attachment nodes
    elemsLocalIDs = [nodesLocalIDs[i] == 1 ?
                     nodesLocalIDs[i] :
                     nodesLocalIDs[i]-1 for i in 1:N]

    # Array of sides of the spring's attachment nodes on the beam elements
    nodesSides = [nodesLocalIDs[i] == 1 ? 
                  1 :
                  2 for i in 1:N]
    
    # Stiffness vectors, resolved in basis s (local basis of the spring)
    kuVec = SVector(ku, 0, 0)
    kpVec = SVector(kp, 0, 0)

    # Pre-deformation vectors, resolved in basis s
    Δu0Vec = SVector(Δu0, 0, 0)
    Δp0Vec = SVector(Δp0, 0, 0)

    # Rotation tensor from basis A to basis b
    if basis == :b
        R0 = if hasDoubleAttachment
            beams[1].geometry.R0
        else
            nodesSides[1] == 1 ? 
            beams[1].elements[elemsLocalIDs[1]].geometry.R0_n1 : 
            beams[1].elements[elemsLocalIDs[1]].geometry.R0_n2
        end
    else#if basis == :A
        R0 = I3
    end

    # Rotation tensor from basis b to basis s
    Rs = rotation_tensor_E321(p)

    # Set stiffness matrices, resolved in basis A
    Ku = SMatrix{3,3,Float64,9}((R0*Rs) * Diagonal(kuVec) * (R0*Rs)')
    Kp = SMatrix{3,3,Float64,9}((R0*Rs) * Diagonal(kpVec) * (R0*Rs)')

    # Pre-deformation vectors, resolved in basis A
    Δu0VecA = SVector{3,Float64}((R0*Rs) * Δu0Vec)
    Δp0VecA = SVector{3,Float64}((R0*Rs) * Δp0Vec)

    # Set force and moment vectors as functions of deformation vectors, and their jacobians, resolved in basis A
    FsFun = SpringForceFun(Δu -> forceFunction(Ku, Δu - Δu0VecA))
    MsFun = SpringForceFun(Δp -> momentFunction(Kp, rotation_between_WM(-Δp0VecA, Δp - Δp0VecA)))

    Fs_ΔuFun = SpringJacFun(Δu -> jacobian3x3_fd(FsFun, SVector{3,Float64}(Δu)))
    Ms_ΔpFun = SpringJacFun(Δp -> jacobian3x3_fd(MsFun, SVector{3,Float64}(Δp)))

    # Set beam construction IDs
    beamConstructionIDs = [beam.constructionID for beam in beams]

    # Update spring assembly data
    setfield!(spring, :assemblyData, SpringAssemblyData(spring,hasDoubleAttachment,beamConstructionIDs,nodesLocalIDs,elemsLocalIDs,nodesSides))
    setfield!(spring, :assemblyConstitutive, SpringAssembledConstitutiveData(Ku,Kp,FsFun,MsFun,Fs_ΔuFun,Ms_ΔpFun))

end