#

#     HingeAxisConstraint composite type

#
@with_kw mutable struct HingeAxisConstraint

    # Primary (inputs to hinge axis constraint creation)
    beam::Beam
    masterElementLocalID::Int64
    slaveElementLocalID::Int64
    localHingeAxis::Vector{<:Real}
    updateHingeAxis::Bool
    loadBalanceLocalNode::Int64
    foldGuessValue::Union{Real,Nothing}
    ΔpValue::Union{Real,Nothing}

    # Secondary (outputs from hinge axis constraint creation)
    initialHingeAxis::Vector{Float64}
    currentHingeAxis::Vector{Float64}
    masterDOF::Int64
    slaveDOFs::Vector{Int64}
    masterElementGlobalID::Int64
    slaveElementGlobalID::Int64
    masterElementGlobalMasterDOF::Int64
    masterElementGlobalSlaveDOFs::Vector{Int64}
    slaveElementGlobalMasterDOF::Int64
    slaveElementGlobalSlaveDOFs::Vector{Int64}
    balanceMoment::Vector{Float64}
    balanceMomentTypes::Vector{String}
    Δp::Vector{Float64}
    Δϕ::Real
    masterElementGlobalDOFs::Vector{Int64}
    slaveElementGlobalDOFs::Vector{Int64}

end
export HingeAxisConstraint


"""
    HingeAxisConstraint(; kwargs...)

Hinge axis constraint constructor

# Keyword arguments
- `beam::Beam` = beam
- `masterElementLocalID::Int64` = local ID of the master (reference) element
- `slaveElementLocalID::Int64` = local ID of the slave element
- `localHingeAxis::Vector{<:Real}` = three-element vector defining the hinge axis, resolved in the undeformed (b) basis of the beam
- `updateHingeAxis::Bool` = flag to update the hinge axis (upon beam deformation)
- `loadBalanceLocalNode::Int64` = local node of the beam where to set the balance load (i.e., the moment necessary to enforce the constraint)
- `foldGuessValue::Union{Real,Nothing}` = initial guess value for the difference in rotation parameter of the master DOF between slave and master elements
- `ΔpValue::Union{Real,Nothing}` = constrained value for the value of the difference between rotation parameters vectors of slave and master elements (value of the rotation across the hinge)
"""
function create_HingeAxisConstraint(; beam::Beam,masterElementLocalID::Int64,slaveElementLocalID::Int64,localHingeAxis::Vector{<:Real},updateHingeAxis::Bool=false,loadBalanceLocalNode::Int64,foldGuessValue::Union{Real,Nothing}=nothing,ΔpValue::Union{Real,Nothing}=nothing)

    # Validate
    @assert masterElementLocalID <= beam.nElements "beam does not have 'masterElementLocalID' elements"
    @assert slaveElementLocalID <= beam.nElements "beam does not have 'slaveElementLocalID' elements"
    @assert length(localHingeAxis) == 3 "specify localHingeAxis as a three-element vector"
    @assert norm(localHingeAxis) > 0 "localHingeAxis can't be a null vector"

    # Normalize local hinge axis
    localHingeAxis /= norm(localHingeAxis)

    # Set initial (undeformed beam) and current hinge axes resolved in basis A
    currentHingeAxis = initialHingeAxis = beam.R0*localHingeAxis

    # Set local master (free) DOF as the greatest direction of local hinge axis
    masterDOF = argmax(abs.(localHingeAxis))

    # Set local slave DOFs
    slaveDOFs = [1,2,3]
    popat!(slaveDOFs,masterDOF)

    # Initialize global IDs (updated later on model assembly)
    masterElementGlobalID, slaveElementGlobalID = 0, 0

    # Initialize global DOFs (updated later on model assembly)
    masterElementGlobalMasterDOF, slaveElementGlobalMasterDOF = 0, 0
    masterElementGlobalSlaveDOFs, slaveElementGlobalSlaveDOFs = [0, 0], [0, 0]
    masterElementGlobalDOFs, slaveElementGlobalDOFs = zeros(Int64,3), zeros(Int64,3)

    # Initialize the balance load vector (moment necessary to enforce the constraint) and the corresponding types (directions)
    balanceMoment = isnothing(ΔpValue) ? zeros(2) : zeros(3)
    if masterDOF == 1
        balanceMomentTypes = isnothing(ΔpValue) ? ["M2A","M3A"] : ["M2A","M3A","M1A"]
    elseif masterDOF == 2
        balanceMomentTypes = isnothing(ΔpValue) ? ["M1A","M3A"] : ["M1A","M3A","M2A"]
    elseif masterDOF == 3
        balanceMomentTypes = isnothing(ΔpValue) ? ["M1A","M2A"] : ["M1A","M2A","M3A"]
    end

    # Initialize vector of rotation parameters across the hinge, resolved in basis A, and angle of rotation [rad]
    Δp, Δϕ = zeros(3), 0.0

    return HingeAxisConstraint(beam=beam,masterElementLocalID=masterElementLocalID,slaveElementLocalID=slaveElementLocalID,localHingeAxis=localHingeAxis,updateHingeAxis=updateHingeAxis,loadBalanceLocalNode=loadBalanceLocalNode,foldGuessValue=foldGuessValue,ΔpValue=ΔpValue,initialHingeAxis=initialHingeAxis,currentHingeAxis=currentHingeAxis,masterDOF=masterDOF,slaveDOFs=slaveDOFs,masterElementGlobalID=masterElementGlobalID, slaveElementGlobalID=slaveElementGlobalID,masterElementGlobalMasterDOF=masterElementGlobalMasterDOF,masterElementGlobalSlaveDOFs=masterElementGlobalSlaveDOFs,slaveElementGlobalMasterDOF=slaveElementGlobalMasterDOF,slaveElementGlobalSlaveDOFs=slaveElementGlobalSlaveDOFs,balanceMoment=balanceMoment,balanceMomentTypes=balanceMomentTypes,Δp=Δp,Δϕ=Δϕ,masterElementGlobalDOFs=masterElementGlobalDOFs,slaveElementGlobalDOFs=slaveElementGlobalDOFs)

end
export create_HingeAxisConstraint