## HINGE AXIS CONSTRAINT COMPOSITE TYPES
# ------------------------------------------------------------------------------
# HingeAxisConstraint assembly data composite type
# Fields
# - `masterGlobalID`: global ID of the master element
# - `slaveGlobalID`: global ID of the slave element
# - `DOF_pM`: global indices of rotation parameters of the master element
# - `DOF_pS`: global indices of rotation parameters of the slave element
# - `eqs_Fp`: global indices of the moment equilibrium equation of the hinge node
# - `DOF_pMS`: global indices of rotation parameters of the master and slave elements
# - `eqs_λ`: global indices of the constraint equations (Lagrange multipliers)
# - `DOF_λ`: global indices of the Lagrange multipliers
#
struct HingeAxisConstraintAssemblyData{N}

    masterGlobalID::Int
    slaveGlobalID::Int

    DOF_pM::SVector{3,Int}
    DOF_pS::SVector{3,Int}

    eqs_Fp::SVector{6,Int}
    DOF_pMS::SVector{6,Int}

    eqs_λ::SVector{N,Int}
    DOF_λ::SVector{N,Int}

end
HingeAxisConstraintAssemblyData{N}() where {N} = 
    HingeAxisConstraintAssemblyData(
        ntuple(_ -> 0, 2)...,
        ntuple(_ -> zeros(SVector{3,Int}), 2)..., 
        ntuple(_ -> zeros(SVector{6,Int}), 2)...,
        ntuple(_ -> zeros(SVector{N,Int}), 2)...,
    )
HingeAxisConstraintAssemblyData(c, masterGlobalID::Int,slaveGlobalID::Int,eqs_Fp::AbstractVector{Int},DOF_pM::AbstractVector{Int},DOF_pS::AbstractVector{Int},eqs_λ::AbstractVector{Int},DOF_λ::AbstractVector{Int}) = begin

    N = length(c.slaveDir)
    DOF_pM = SVector{3,Int}(DOF_pM)
    DOF_pS = SVector{3,Int}(DOF_pS)
    eqs_Fp = SVector{6,Int}(eqs_Fp)
    DOF_pMS = SVector{6,Int}(vcat(DOF_pM,DOF_pS))
    eqs_λ = SVector{N,Int}(eqs_λ)
    DOF_λ = SVector{N,Int}(DOF_λ)

    HingeAxisConstraintAssemblyData{N}(
        masterGlobalID,
        slaveGlobalID,
        DOF_pM,
        DOF_pS,
        eqs_Fp,
        DOF_pMS,
        eqs_λ,
        DOF_λ,
    )

end

# HingeAxisConstraint Jacobian indices maps composite type
# Fields
# - `Fp_pMS`: moment equilibrium equations and rotation parameters of the master and slave elements
# - `Fp_λ`: moment equilibrium equations and Lagrange multipliers
# - `Fλ_pMS`: constraint equations and rotation parameters of the master and slave elements
struct HingeAxisConstraintJacobianIndicesMaps{N,N6}

    Fp_pMS::SVector{36,Int}
    Fp_λ::SVector{N6,Int}
    Fλ_pMS::SVector{N6,Int}

end
HingeAxisConstraintJacobianIndicesMaps{N,N6}() where {N,N6} = 
    HingeAxisConstraintJacobianIndicesMaps{N,6*N}(
        zeros(SVector{36,Int}),
        zeros(SVector{6*N,Int}),
        zeros(SVector{6*N,Int}),
    )

# HingeAxisConstraint state composite type
# Fields
# - `λ`: Lagrange multipliers vector
# - `B`: ∂C/∂p, where C is the constraint equations' vector and p is the vector of rotation parameters of the master and slave elements (constraint matrix or Jacobian of the constraints)
# - `Bᵀ` = transpose of B
# - `Gᵀ`: mapping from virtual rotations to variations of the rotation parameters of the master and slave elements (combined G matrix)
# - `GᵀBᵀ`: (B*G)ᵀ, where G maps virtual rotations (resolved in basis A) into variations of the rotation parameters of the master and slave elements (δp = G δθ); GᵀBᵀ*λ are the constraint moments, resolved in basis A, acting on the master (first three entries) and slave (last three entries) elements
# - `GᵀBᵀλ_p`: ∂(GᵀBᵀ*λ)/∂p
# - `pM`: master element rotation parameters vector
# - `pS`: slave element rotation parameters vector
# - `pH`: hinge rotation parameters vector
# - `deformedHingeAxis`: deformed hinge axis resolved in basis A
# - `balanceMoment`: balance moment vector resolved in basis A
# - `pHValue`: constrained value for the hinge rotation
# - `ϕ`: hinge rotation angle
# - `hingeMoment`: hinge moment magnitude
#
mutable struct HingeAxisConstraintState{N}

    λ::MVector{N,Float64}

    B::MMatrix{N,6,Float64}
    Bᵀ::MMatrix{6,N,Float64}
    Gᵀ::MMatrix{6,6,Float64}

    GᵀBᵀ::MMatrix{6,N,Float64}
    GᵀBᵀλ_p::MMatrix{6,6,Float64}
    
    pM::MVector{3,Float64}
    pS::MVector{3,Float64}
    pH::MVector{3,Float64}

    deformedHingeAxis::MVector{3,Float64}
    balanceMoment::MVector{3,Float64}

    pHValue::Float64
    ϕ::Float64
    hingeMoment::Float64

end
HingeAxisConstraintState{N}() where {N} = 
    HingeAxisConstraintState(
        zeros(MVector{N,Float64}),
        zeros(MMatrix{N,6,Float64}), zeros(MMatrix{6,N,Float64}), zeros(MMatrix{6,6,Float64}),
        zeros(MMatrix{6, N, Float64}), zeros(MMatrix{6,6,Float64}),
        ntuple(_ -> zeros(MVector{3,Float64}), 3)...,
        zeros(MVector{3,Float64}), zeros(MVector{3,Float64}),
        ntuple(_ -> 0.0, 3)...,
    )
Base.copy(x::HingeAxisConstraintState{N}) where {N} = copy_fields(x)
const HACStateUnion = Union{HingeAxisConstraintState{2}, HingeAxisConstraintState{3}}

# HingeAxisConstraint composite type
# Fields
# - `beam`: beam to which the hinge axis constraint is applied
# - `masterOnLeft`: flag for the master element being on the left of the slave
# - `localHingeAxis`: hinge axis, resolved in the undeformed (b) basis of the beam
# - `pHValuePrescribed`: constrained value for the rotation at the hinge
# - `rotationIsFixed`: flag for hinge rotation magnitude being known (fixed)
# - `masterLocalID`: local ID of the master element
# - `slaveLocalID`: local ID of the slave element
# - `masterDir`: direction of the master element's rotation DOF
# - `slaveDir`: directions of the slave element's rotation DOFs
# - `initialHingeAxis`: initial (undeformed beam) hinge axis resolved in basis A
# - `assemblyData`: assembly data for the hinge axis constraint
# - `state`: state of the hinge axis constraint
#
@with_kw mutable struct HingeAxisConstraint{N,N6}

    # Definition fields (inputs for HingeAxisConstraint creation)
    beam::Beam
    masterOnLeft::Bool
    localHingeAxis::SVector{3,Float64}
    pHValuePrescribed::Union{Missing,Float64}

    # Deduced fields (from HingeAxisConstraint inputs)
    rotationIsFixed::Bool
    masterLocalID::Int
    slaveLocalID::Int
    masterDir::Int
    slaveDir::SVector{N,Int}
    initialHingeAxis::SVector{3,Float64}

    # Assembly-dependent fields (to be set by Model)
    assemblyData::HingeAxisConstraintAssemblyData{N} = HingeAxisConstraintAssemblyData{N}()

    # Jacobian indices maps for the hinge axis constraint
    jacobianMaps::HingeAxisConstraintJacobianIndicesMaps{N,N6} = HingeAxisConstraintJacobianIndicesMaps{N,N6}()

    # State field
    state::HingeAxisConstraintState{N} = HingeAxisConstraintState{N}()
    
end
export HingeAxisConstraint


## PUBLIC API
# ------------------------------------------------------------------------------
"""
    create_HingeAxisConstraint(; kwargs...)

Hinge axis constraint constructor

# Keyword arguments
- `beam::Beam`: beam
- `masterOnLeft::Bool`: flag for the master element being on the left of the slave element
- `localHingeAxis::AbstractVector{<:Real}`: three-element vector defining the hinge axis, resolved in the undeformed (b) basis of the beam
- `pHValue::Union{Nothing,Real}`: constrained value for the Wiener-Milenkovic rotation at the hinge
"""
function create_HingeAxisConstraint(; 
    beam::Beam,
    masterOnLeft::Bool=true,
    localHingeAxis::AbstractVector{<:Real},
    pHValue::Union{Missing,Real}=missing,
    )

    # Validate
    @assert (!isempty(beam.hinges.hingedNodes) && length(beam.hinges.hingedNodes) == 1) "beam must have only one hinged node"
    @assert length(localHingeAxis) == 3 "specify localHingeAxis as a three-element vector"
    @assert norm(localHingeAxis) > 0 "localHingeAxis can't be a null vector"

    # Flag for hinge rotation magnitude being known (fixed)
    rotationIsFixed = !ismissing(pHValue)

    # Number of unknowns in the constraint
    N = rotationIsFixed ? 
        3 :
        2

    # Set local IDs of master and slave elements
    hingeNode = beam.hinges.hingedNodes[1]
    masterLocalID = masterOnLeft ? 
                    hingeNode-1 :
                    hingeNode
    slaveLocalID = masterOnLeft ?
                   hingeNode :
                   hingeNode-1

    # Normalize local hinge axis
    localHingeAxis /= norm(localHingeAxis)

    # Set initial (undeformed beam) hinge axis resolved in basis A
    initialHingeAxis = beam.geometry.R0 * localHingeAxis

    # Set master (free) and slave directions
    slaveDir = [1,2,3]
    if rotationIsFixed
        masterDir = 0
    else
        masterDir = argmax(abs.(initialHingeAxis))
        popat!(slaveDir,masterDir)
    end

    # Change variable name and set type for pHValue 
    pHValuePrescribed = ismissing(pHValue) ? 
                        missing : 
                        Float64(pHValue)

    # Set correct types
    localHingeAxis = SVector{3,Float64}(localHingeAxis)
    initialHingeAxis = SVector{3,Float64}(initialHingeAxis)
    slaveDir = SVector{N,Int}(slaveDir)

    return HingeAxisConstraint{N,6*N}(
        beam=beam,
        masterOnLeft=masterOnLeft,
        localHingeAxis=localHingeAxis,
        pHValuePrescribed=pHValuePrescribed,
        rotationIsFixed=rotationIsFixed,
        masterLocalID=masterLocalID,
        slaveLocalID=slaveLocalID,
        masterDir=masterDir,
        slaveDir=slaveDir,
        initialHingeAxis=initialHingeAxis,
        )

end
export create_HingeAxisConstraint


## PRIVATE API
# ------------------------------------------------------------------------------
# Hinge constraint equations
function C(pM,pS,initialHingeAxis,rotationIsFixed,pHValuePrescribed,slaveDir)

    # Scale rotation parameters of the slave element when rotation is not fixed(pHValue is not prescribed)
    pS = rotationIsFixed ? 
         pS : 
         scaled_rotation_parameters(pS)

    # Hinge rotation parameters
    pH = hinge_rotation_parameters(pM,pS)

    # Deformed hinge axis
    n = deformed_hinge_axis(pM,initialHingeAxis)

    # Hinge rotation value
    ψ = coalesce(pHValuePrescribed, hinge_rotation_value(pH,n))

    # Constraint equation
    c = pH - ψ * n

    return c[slaveDir]

end

# Computes the hinge rotation parameters, given the master and slave rotation parameters
function hinge_rotation_parameters(pM,pS)

    # Master and slave rotation tensors
    RM = rotation_tensor_WM(pM)
    RS = rotation_tensor_WM(pS)

    # Hinge rotation tensor
    RH = RS * RM'

    # Hinge rotation parameters
    pH = rotation_parameters_WM(RH)

    return pH

end

# Computes the deformed hinge axis, given the master element rotation parameters (pM) and the initial hinge axis (n₀)
function deformed_hinge_axis(pM,n₀)
    return rotation_tensor_WM(pM) * n₀
end

# Computes the signed magnitude of the rotation, given the hinge rotation parameters (pH) and the deformed hinge axis (n)
function hinge_rotation_value(pH,n)
    return dot(pH, n)
end

# Computes the slave rotation parameters, given the master rotation parameters (pM), the initial hinge axis (n₀), and the rotation magnitude (pHValue)
function pS_from_pM_and_pHValue(pM,n₀,pHValue)

    # Current hinge rotation axis
    n = deformed_hinge_axis(pM,n₀)

    # Hinge rotation parameters and associated rotation tensor
    pH = pHValue * n
    RH = rotation_tensor_WM(pH)

    # Master element rotation tensor
    RM = rotation_tensor_WM(pM)

    # Slave element rotation tensor and associated parameters
    RS = RH * RM
    pS = rotation_parameters_WM(RS)

    return pS

end

# Jacobian of constraint equations w.r.t. system states (constraint matrix B)
function ∂C_∂p(pM,pS,initialHingeAxis,rotationIsFixed,pHValuePrescribed,slaveDir)
    return ForwardDiff.jacobian(x -> C(SVector{3}(x[1],x[2],x[3]), SVector{3}(x[4],x[5],x[6]),initialHingeAxis,rotationIsFixed,pHValuePrescribed,slaveDir), SVector{6}(pM[1],pM[2],pM[3],pS[1],pS[2],pS[3]))
end

# Computes the transpose of the compound matrix G = [GM 0; 0 GS], where GM and GS are the matrices for the master and slave elements that map virtual rotations (resolved in basis A) into variations of the (extended) rotation parameters: δpM = GM δθM and δpS = GS δθS
function GT(pM,pS)

    # Convert to static vectors
    pM = SVector{3}(pM)
    pS = SVector{3}(pS)

    # Compute the G matrices for the master and slave elements
    GM = HTinv_WM(pM) * rotation_tensor_WM(pM)'
    GS = HTinv_WM(pS) * rotation_tensor_WM(pS)'

    # Gᵀ = [GMᵀ 0; 0 GSᵀ]
    return vcat(hcat(GM',Z3), hcat(Z3,GS'))

end

# Derivative of the constraint moments ( (B*G)ᵀ*λ ) w.r.t. the master and slave rotation parameters
function ∂GTBTλ_∂p(λ,Bᵀ,Gᵀ,pM,pS,initialHingeAxis,rotationIsFixed,pHValuePrescribed,slaveDir)

    # State vector containing the master and slave rotation parameters
    x = SVector{6}(pM[1],pM[2],pM[3],pS[1],pS[2],pS[3])

    # Constraint forces conjugate to the rotation parameters
    Bᵀλ = SVector{6}(Bᵀ * λ)

    # ∂Gᵀ/∂p * Bᵀλ
    Gᵀv_p = ForwardDiff.jacobian(x -> GT(SVector{3}(x[1],x[2],x[3]), SVector{3}(x[4],x[5],x[6])) * Bᵀλ, x)

    # ∂(Bᵀλ)/∂p = ∇²(λᵀC)
    Bᵀλ_p = ForwardDiff.hessian(x -> dot(λ, C(SVector{3}(x[1],x[2],x[3]), SVector{3}(x[4],x[5],x[6]), initialHingeAxis,rotationIsFixed,pHValuePrescribed,slaveDir)), x)

    return Gᵀv_p  + Gᵀ * Bᵀλ_p

end