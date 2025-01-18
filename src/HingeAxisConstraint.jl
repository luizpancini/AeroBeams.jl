#

#     HingeAxisConstraint composite type

#
@with_kw mutable struct HingeAxisConstraint

    # Primary (inputs to hinge axis constraint creation)
    beam::Beam
    masterElementLocalID::Int64
    slaveElementLocalID::Int64
    localHingeAxis::Vector{<:Real}
    loadBalanceLocalNode::Int64
    pHValue::Union{Real,Nothing}
    pHGuessValue::Union{Real,Nothing}

    # Secondary (outputs from hinge axis constraint creation)
    rotationIsFixed::Bool
    initialHingeAxis::Vector{Float64}
    currentHingeAxis::Vector{Float64}
    masterDOF::Int64
    slaveDOFs::Vector{Int64}
    masterElementGlobalID::Int64
    slaveElementGlobalID::Int64
    masterElementGlobalDOFs::Vector{Int64}
    slaveElementGlobalDOFs::Vector{Int64}
    balanceMoment::Vector{Float64}
    balanceMomentDirections::Vector{String}
    pH::Vector{Float64}
    ϕ::Real

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
- `loadBalanceLocalNode::Int64` = local node of the beam where to set the balance load (i.e., the moment necessary to enforce the constraint)
- `pHValue::Union{Real,Nothing}` = constrained value for the rotation at the hinge
- `pHGuessValue::Union{Real,Nothing}` = initial guess value for the rotation at the hinge (in case it is unknonwn)
"""
function create_HingeAxisConstraint(; beam::Beam,masterElementLocalID::Int64,slaveElementLocalID::Int64,localHingeAxis::Vector{<:Real},loadBalanceLocalNode::Int64,pHValue::Union{Real,Nothing}=nothing,pHGuessValue::Union{Real,Nothing}=nothing)

    # Validate
    @assert masterElementLocalID <= beam.nElements "beam does not have 'masterElementLocalID' elements"
    @assert slaveElementLocalID <= beam.nElements "beam does not have 'slaveElementLocalID' elements"
    @assert length(localHingeAxis) == 3 "specify localHingeAxis as a three-element vector"
    @assert norm(localHingeAxis) > 0 "localHingeAxis can't be a null vector"
    @assert !(!isnothing(pHValue) && !isnothing(pHGuessValue)) "if pHValue is input, do not input pHGuessValue (and vice-versa)"

    # TF for hinge rotation magnitude being known (fixed)
    rotationIsFixed = !isnothing(pHValue)

    # Normalize local hinge axis
    localHingeAxis /= norm(localHingeAxis)

    # Set initial (undeformed beam) and current hinge axes resolved in basis A
    currentHingeAxis = initialHingeAxis = beam.R0*localHingeAxis

    # Set master (free) and slave DOFs in case of unknonwn hinge rotation
    if rotationIsFixed
        masterDOF = 0
        slaveDOFs = zeros(Int64,2)
    else
        masterDOF = argmax(abs.(initialHingeAxis))
        slaveDOFs = [1,2,3]
        popat!(slaveDOFs,masterDOF)
    end

    # Initialize global IDs (updated later on model assembly)
    masterElementGlobalID, slaveElementGlobalID = 0, 0

    # Initialize global DOFs (updated later on model assembly)
    masterElementGlobalDOFs, slaveElementGlobalDOFs = zeros(Int64,3), zeros(Int64,3)

    # Initialize the balance moment vector (moment necessary to enforce the constraint) and corresponding directions
    balanceMoment = rotationIsFixed ? zeros(3) : zeros(2)
    balanceMomentDirections = rotationIsFixed ? ["M1A","M2A","M3A"] : ifelse(masterDOF==1,["M2A","M3A"],ifelse(masterDOF==2,["M1A","M3A"],["M1A","M2A"]))

    # Initialize vector of rotation parameters across the hinge, resolved in basis A, and angle of rotation [rad]
    pH, ϕ = zeros(3), 0.0

    return HingeAxisConstraint(beam=beam,masterElementLocalID=masterElementLocalID,slaveElementLocalID=slaveElementLocalID,localHingeAxis=localHingeAxis,loadBalanceLocalNode=loadBalanceLocalNode,pHValue=pHValue,pHGuessValue=pHGuessValue,rotationIsFixed=rotationIsFixed,initialHingeAxis=initialHingeAxis,currentHingeAxis=currentHingeAxis,masterDOF=masterDOF,slaveDOFs=slaveDOFs,masterElementGlobalID=masterElementGlobalID, slaveElementGlobalID=slaveElementGlobalID,masterElementGlobalDOFs=masterElementGlobalDOFs,slaveElementGlobalDOFs=slaveElementGlobalDOFs,balanceMoment=balanceMoment,balanceMomentDirections=balanceMomentDirections,pH=pH,ϕ=ϕ)

end
export create_HingeAxisConstraint