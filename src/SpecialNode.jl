## SPECIAL NODE COMPOSITE TYPES
# ------------------------------------------------------------------------------
# SpecialNode boundary condition flags composite type
# Fields
# - `hasFollowers`: boolean indicating if the special node has follower forces
# - `allBCflags`: vector of tuples containing the isLoad, isTrim and isFollower flags accounting for all BCs applied to the special node
# - `uIsPrescribed`: boolean vector for translational displacements being prescribed
# - `pIsPrescribed`: boolean vector for rotational displacements being prescribed
# - `trimIsPrescribed`: boolean vector for trim loads being prescribed
#
struct BCFlags

    hasFollowers::Bool

    allBCflags::Vector{Tuple{SVector{6,Bool},SVector{6,Bool},SVector{6,Bool}}}

    uIsPrescribed::SVector{3,Bool}
    pIsPrescribed::SVector{3,Bool}
    trimIsPrescribed::SVector{6,Bool}

    uNotPrescribedMat::SMatrix{3,3,Bool,9}
    pNotPrescribedMat::SMatrix{3,3,Bool,9}

end


# SpecialNode assembly data composite type
# Fields
# - `globalID`: global ID of the special node
# - `connectedElementsGlobalIDs`: global IDs of the connected elements
# - `ζonElements`: ζ coordinates of the special node on the connected elements
# - `doubleAttachmentSpringsData`: data for double attachment springs connected to the special node
# - `eqs_Fu`: global indices of the force equilibrium equations
# - `eqs_Fp`: global indices of the moment equilibrium equations
# - `eqs_FF`: global indices of the translational displacements compatibility equations (for each connected element)
# - `eqs_FM`: global indices of the rotational displacements compatibility equations (for each connected element)
# - `DOF_uF`: global DOFs of the translational displacements and forces
# - `DOF_pM`: global DOFs of the rotational displacements and moments
# - `DOF_T`: global DOFs of the trim loads
#
struct SpecialNodeAssemblyData

    globalID::Int
    connectedElementsGlobalIDs::Vector{Int}
    ζonElements::Vector{Int}

    doubleAttachmentSpringsData::Vector{DoubleAttachmentSpringData}

    eqs_Fu::SVector{3,Int}
    eqs_Fp::SVector{3,Int}
    eqs_FF::Vector{SVector{3,Int}}
    eqs_FM::Vector{SVector{3,Int}}

    DOF_uF::SVector{3,Int}
    DOF_pM::SVector{3,Int}
    DOF_T::SVector{6,Int}

end
SpecialNodeAssemblyData(globalID,connectedElementsGlobalIDs,ζonElements) =
    SpecialNodeAssemblyData(
        globalID, connectedElementsGlobalIDs, ζonElements,
        Vector{DoubleAttachmentSpringData}(undef, 0),
        zeros(SVector{3,Int}), zeros(SVector{3,Int}), Vector{SVector{3,Int}}(undef, 0), Vector{SVector{3,Int}}(undef, 0),
        zeros(SVector{3,Int}), zeros(SVector{3,Int}), zeros(SVector{6,Int}),
    )

SpecialNodeAssemblyData(specialNode, DOF_uF::SVector{3,Int},DOF_pM::SVector{3,Int}) = 
    SpecialNodeAssemblyData(
        specialNode.assemblyData.globalID,
        specialNode.assemblyData.connectedElementsGlobalIDs,
        specialNode.assemblyData.ζonElements,
        specialNode.assemblyData.doubleAttachmentSpringsData,
        specialNode.assemblyData.eqs_Fu,
        specialNode.assemblyData.eqs_Fp,
        specialNode.assemblyData.eqs_FF,
        specialNode.assemblyData.eqs_FM,
        DOF_uF,
        DOF_pM,
        specialNode.assemblyData.DOF_T
    )

SpecialNodeAssemblyData(specialNode, DOF_T::SVector{6,Int}) = 
    SpecialNodeAssemblyData(
        specialNode.assemblyData.globalID,
        specialNode.assemblyData.connectedElementsGlobalIDs,
        specialNode.assemblyData.ζonElements,
        specialNode.assemblyData.doubleAttachmentSpringsData,
        specialNode.assemblyData.eqs_Fu,
        specialNode.assemblyData.eqs_Fp,
        specialNode.assemblyData.eqs_FF,
        specialNode.assemblyData.eqs_FM,
        specialNode.assemblyData.DOF_uF,
        specialNode.assemblyData.DOF_pM,
        DOF_T
    )

SpecialNodeAssemblyData(specialNode, eqs_Fu::SVector{3,Int},eqs_Fp::SVector{3,Int},eqs_FF::AbstractVector{SVector{3,Int}},eqs_FM::AbstractVector{SVector{3,Int}}) = 
    SpecialNodeAssemblyData(
        specialNode.assemblyData.globalID,
        specialNode.assemblyData.connectedElementsGlobalIDs,
        specialNode.assemblyData.ζonElements,
        specialNode.assemblyData.doubleAttachmentSpringsData,
        eqs_Fu,
        eqs_Fp,
        eqs_FF,
        eqs_FM,
        specialNode.assemblyData.DOF_uF,
        specialNode.assemblyData.DOF_pM,
        specialNode.assemblyData.DOF_T
    )

SpecialNodeAssemblyData(specialNode, doubleAttachmentSpringsData::Vector{DoubleAttachmentSpringData}) =
    SpecialNodeAssemblyData(
        specialNode.assemblyData.globalID,
        specialNode.assemblyData.connectedElementsGlobalIDs,
        specialNode.assemblyData.ζonElements,
        doubleAttachmentSpringsData,
        specialNode.assemblyData.eqs_Fu,
        specialNode.assemblyData.eqs_Fp,
        specialNode.assemblyData.eqs_FF,
        specialNode.assemblyData.eqs_FM,
        specialNode.assemblyData.DOF_uF,
        specialNode.assemblyData.DOF_pM,
        specialNode.assemblyData.DOF_T
    )


# SpecialNode Jacobian indices maps composite type
# Fields
# - `Fu_uF`: indices for equilibrium equations corresponding to translational DOFs
# - `Fp_pM`: indices for equilibrium equations corresponding to rotational DOFs
# - `Fu_uFother`: indices for equilibrium equations corresponding to translational DOFs of other spring attachment nodes
# - `Fp_pMother`: indices for equilibrium equations corresponding to rotational DOFs of other spring attachment nodes
#
struct SpecialNodeJacobianIndicesMaps

    Fu_uF::SVector{9,Int}
    Fp_pM::SVector{9,Int}

    Fu_uFother::Vector{SVector{9,Int}}
    Fp_pMother::Vector{SVector{9,Int}}

end
SpecialNodeJacobianIndicesMaps() = SpecialNodeJacobianIndicesMaps(
    ntuple(_ -> zeros(SVector{9,Int}), 2)...,
    ntuple(_ -> Vector{SVector{9,Int}}(), 2)...,
)


# SpecialNode state composite type
# Fields
# - `u`: translational displacement vector resolved in basis A
# - `p`: rotation parameters vector from basis b to basis B, resolved in basis A
# - `F`: sectional force vector resolved in basis A
# - `M`: sectional moment vector resolved in basis A
# - `R`: rotation tensor associated with the rotation parameters `p`
# - `F_p`: ∂F/∂p
# - `M_p`: ∂M/∂p
#
struct SpecialNodeState

    u::MVector{3,Float64}
    p::MVector{3,Float64}
    F::MVector{3,Float64}
    M::MVector{3,Float64}

    R::MMatrix{3,3,Float64,9}
    
    F_p::MMatrix{3,3,Float64,9}
    M_p::MMatrix{3,3,Float64,9}

end
SpecialNodeState() = SpecialNodeState(
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
    zeros(MMatrix{3,3,Float64,9}),
    ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 2)...,
)
Base.copy(x::SpecialNodeState) = copy_fields(x)


# SpecialNode composite type
# Fields
# - `BCs`: boundary conditions applied to the special node
# - `springs`: springs attached to the special node
# - `flags`: flags for generalized displacements and trim loads being prescribed
# - `assemblyData`: assembly data for the special node
# - `state`: state of the special node
#
@with_kw mutable struct SpecialNode

    # BCs, springs and flags
    BCs::Vector{BC}
    springs::Vector{Spring}
    flags::BCFlags

    # Assembly-dependent fields (to be set by Model)
    assemblyData::SpecialNodeAssemblyData
    jacobianMaps::SpecialNodeJacobianIndicesMaps = SpecialNodeJacobianIndicesMaps()

    # State field
    state::SpecialNodeState = SpecialNodeState()

end


## PRIVATE API
# ------------------------------------------------------------------------------
# Special node constructor
function create_SpecialNode(globalID,connectedElementsGlobalIDs,ζonElements,BCs,springs)

    # Flags for generalized displacements and trim loads being prescribed 
    flags = special_node_bc_flags(BCs)

    # Initialize assembly data
    assemblyData = SpecialNodeAssemblyData(globalID,connectedElementsGlobalIDs,ζonElements)

    return SpecialNode(
        BCs=BCs,
        springs=springs,
        flags=flags,
        assemblyData=assemblyData,
        )

end


# Sets flags for generalized displacements and trim loads being prescribed based on the BCs applied to the special node
function special_node_bc_flags(BCs)

    # Initialize flags for DOFs being loads and trim loads
    isLoad = trues(6)
    isTrim = falses(6)

    # Loop over BCs and update flags
    for BC in BCs
        isLoad = isLoad .& BC.deducedData.isLoad
        isTrim = isTrim .| BC.deducedData.isTrim
    end

    # Set flag for follower forces
    hasFollowers = any(vcat([BC.deducedData.isFollower for BC in BCs]...))

    # Set flags for DOFs being loads, trim loads and follower loads for all BCs applied to the special node 
    allBCflags = [(BC.deducedData.isLoad, BC.deducedData.isTrim, BC.deducedData.isFollower) for BC in BCs]

    # Set flags for DOFs being prescribed (not loads) and trim loads
    uIsPrescribed = .!isLoad[1:3]
    pIsPrescribed = .!isLoad[4:6]
    trimIsPrescribed = isTrim

    # Create diagonal matrices indicating which DOFs are not prescribed
    uNotPrescribedMat = Diagonal(.!uIsPrescribed)
    pNotPrescribedMat = Diagonal(.!pIsPrescribed)

    # Set correct types
    uIsPrescribed = SVector{3,Bool}(uIsPrescribed)
    pIsPrescribed = SVector{3,Bool}(pIsPrescribed)
    trimIsPrescribed = SVector{6,Bool}(trimIsPrescribed)
    uNotPrescribedMat = SMatrix{3,3,Bool,9}(uNotPrescribedMat)
    pNotPrescribedMat = SMatrix{3,3,Bool,9}(pNotPrescribedMat)

    return BCFlags(
        hasFollowers,
        allBCflags,
        uIsPrescribed,
        pIsPrescribed,
        trimIsPrescribed,
        uNotPrescribedMat,
        pNotPrescribedMat,
    )

end