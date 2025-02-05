#

#     HingeAxisConstraint composite type

#
@with_kw mutable struct HingeAxisConstraint

    # Primary (inputs to hinge axis constraint creation)
    beam::Beam
    localHingeAxis::Vector{<:Real}
    pHValue::Union{Real,Nothing}
    solutionMethod::String

    # Secondary (outputs from hinge axis constraint creation)
    rotationIsFixed::Bool
    initialHingeAxis::Vector{Float64}
    currentHingeAxis::Vector{Float64}
    masterElementLocalID::Int64
    slaveElementLocalID::Int64
    masterDOF::Int64
    slaveDOFs::Vector{Int64}
    masterElementGlobalID::Int64
    slaveElementGlobalID::Int64
    masterElementGlobalDOFs::Vector{Int64}
    slaveElementGlobalDOFs::Vector{Int64}
    balanceMomentNode::Int64
    balanceMoment::Vector{Float64}
    balanceMomentDirections::Vector{String}
    balanceMomentBCID::Int64
    λ::Vector{Float64}
    Jc::Matrix{Float64}
    pH::Vector{Float64}
    ϕ::Real

end
export HingeAxisConstraint


"""
    HingeAxisConstraint(; kwargs...)

Hinge axis constraint constructor

# Keyword arguments
- `beam::Beam`: beam
- `localHingeAxis::Vector{<:Real}`: three-element vector defining the hinge axis, resolved in the undeformed (b) basis of the beam
- `pHValue::Union{Real,Nothing}`: constrained value for the rotation at the hinge
- `solutionMethod::String`: solution method for the constraint ("addedResidual" or "appliedMoment")
"""
function create_HingeAxisConstraint(; beam::Beam,localHingeAxis::Vector{<:Real},pHValue::Union{Real,Nothing}=nothing,solutionMethod::String="appliedMoment")

    # Validate
    @assert solutionMethod in ["addedResidual", "appliedMoment"] "set solutionMethod as 'addedResidual' or 'appliedMoment'"
    @assert (!isempty(beam.hingedNodes) && length(beam.hingedNodes) == 1) "beam must have one hinged node"
    @assert length(localHingeAxis) == 3 "specify localHingeAxis as a three-element vector"
    @assert norm(localHingeAxis) > 0 "localHingeAxis can't be a null vector"

    # Set local IDs of master and slave elements
    hingeNode = beam.hingedNodes[1]
    masterElementLocalID = hingeNode-1
    slaveElementLocalID = hingeNode

    # Set ID of the node where the balance loads will be applied, in case solutionMethod is "appliedMoment"
    balanceMomentNode = hingeNode+1

    # TF for hinge rotation magnitude being known (fixed)
    rotationIsFixed = !isnothing(pHValue)

    # Normalize local hinge axis
    localHingeAxis /= norm(localHingeAxis)

    # Set initial (undeformed beam) and current hinge axes resolved in basis A
    currentHingeAxis = initialHingeAxis = beam.R0*localHingeAxis

    # Set master (free) and slave DOFs in case of unknonwn hinge rotation
    slaveDOFs = [1,2,3]
    if rotationIsFixed
        masterDOF = 0
    else
        masterDOF = argmax(abs.(initialHingeAxis))
        popat!(slaveDOFs,masterDOF)
    end

    # Initialize global IDs (updated later on model assembly)
    masterElementGlobalID, slaveElementGlobalID = 0, 0

    # Initialize global DOFs (updated later on model assembly)
    masterElementGlobalDOFs, slaveElementGlobalDOFs = zeros(Int64,3), zeros(Int64,3)

    # Initialize the balance moment vector (moment necessary to enforce the constraint)
    balanceMoment = rotationIsFixed ? zeros(3) : zeros(2)
    balanceMomentDirections = rotationIsFixed ? ["M1A","M2A","M3A"] : ifelse(masterDOF==1,["M2A","M3A"],ifelse(masterDOF==2,["M1A","M3A"],["M1A","M2A"]))

    # Initialize ID of the BC corresponding to the balance moment
    balanceMomentBCID = 0

    # Initialize the array of Lagrange multipliers
    λ = balanceMoment

    # Initialize Jacobian matrix of the constraint (size is updated later on model assembly)
    Jc = zeros(0,0)

    # Initialize vector of rotation parameters across the hinge, resolved in basis A, and angle of rotation [rad]
    pH, ϕ = zeros(3), 0.0

    return HingeAxisConstraint(beam=beam,localHingeAxis=localHingeAxis,pHValue=pHValue,solutionMethod=solutionMethod,rotationIsFixed=rotationIsFixed,initialHingeAxis=initialHingeAxis,currentHingeAxis=currentHingeAxis,masterElementLocalID=masterElementLocalID,slaveElementLocalID=slaveElementLocalID,masterDOF=masterDOF,slaveDOFs=slaveDOFs,masterElementGlobalID=masterElementGlobalID, slaveElementGlobalID=slaveElementGlobalID,masterElementGlobalDOFs=masterElementGlobalDOFs,slaveElementGlobalDOFs=slaveElementGlobalDOFs,balanceMomentNode=balanceMomentNode,balanceMoment=balanceMoment,balanceMomentDirections=balanceMomentDirections,balanceMomentBCID=balanceMomentBCID,λ=λ,Jc=Jc,pH=pH,ϕ=ϕ)

end
export create_HingeAxisConstraint


#
# mutable struct HingeAxisConstraintData

#     HingeAxisConstraintData composite type

# Fields
# `balanceMoment::Vector{Float64}`: balance moment vector around the hinge, resolved in basis A
# `pH::Vector{Float64}`: hinge rotation parameters, resolved in basis A
# `ϕ::Real`: hinge angle
#
mutable struct HingeAxisConstraintData
    
    # Fields
    balanceMoment::Vector{Float64}
    pH::Vector{Float64}
    ϕ::Real

    # Constructor
    function HingeAxisConstraintData(balanceMoment,pH,ϕ)
        return new(balanceMoment,pH,ϕ)
    end
end