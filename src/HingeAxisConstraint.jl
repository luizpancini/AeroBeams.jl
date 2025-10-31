#

#     HingeAxisConstraint composite type

#
@with_kw mutable struct HingeAxisConstraint

    # Primary (inputs to hinge axis constraint creation)
    beam::Beam
    masterOnLeft::Bool
    localHingeAxis::Vector{<:Real}
    pHValue::Union{Real,Nothing}
    solutionMethod::String
    updateAllDOFinResidual::Bool

    # Secondary (outputs from hinge axis constraint creation)
    rotationIsFixed::Bool
    initialHingeAxis::Vector{Float64}
    deformedHingeAxis::Vector{Float64}
    masterLocalID::Int
    slaveLocalID::Int
    masterDir::Int
    slaveDir::Vector{Int}
    masterGlobalID::Int
    slaveGlobalID::Int
    masterGlobalDOFs::Vector{Int}
    slaveGlobalDOFs::Vector{Int}
    constraintGlobalDOFs::Vector{Int}
    λEqs::Vector{Int}
    balanceMomentNode::Int
    balanceMoment::Vector{Float64}
    balanceMomentDirections::Vector{String}
    balanceMomentBCID::Int
    λ::Vector{Float64}
    Jc::Matrix{Float64}
    pH::Vector{Float64}
    ϕ::Real
    hingeMoment::Real

end
export HingeAxisConstraint


"""
    HingeAxisConstraint(; kwargs...)

Hinge axis constraint constructor

# Keyword arguments
- `beam::Beam`: beam
- `masterOnLeft::Bool`: flag for the master element being on the left of the slave element
- `localHingeAxis::Vector{<:Real}`: three-element vector defining the hinge axis, resolved in the undeformed (b) basis of the beam
- `pHValue::Union{Real,Nothing}`: constrained value for the rotation at the hinge
- `solutionMethod::String`: solution method for the constraint ("addedResidual" or "appliedMoment")
- `updateAllDOFinResidual::Bool`: flag to update all rotation DOFs in the residual array
"""
function create_HingeAxisConstraint(; beam::Beam,masterOnLeft::Bool=true,localHingeAxis::Vector{<:Real},pHValue::Union{Real,Nothing}=nothing,solutionMethod::String="addedResidual",updateAllDOFinResidual::Bool=true)

    # Validate
    @assert solutionMethod in ["addedResidual", "appliedMoment"] "set solutionMethod as 'addedResidual' or 'appliedMoment'"
    @assert (!isempty(beam.hingedNodes) && length(beam.hingedNodes) == 1) "beam must have one hinged node"
    @assert length(localHingeAxis) == 3 "specify localHingeAxis as a three-element vector"
    @assert norm(localHingeAxis) > 0 "localHingeAxis can't be a null vector"

    # Set local IDs of master and slave elements
    hingeNode = beam.hingedNodes[1]
    masterLocalID = masterOnLeft ? hingeNode-1 : hingeNode
    slaveLocalID = masterOnLeft ? hingeNode : hingeNode-1

    # Set ID of the node where the balance loads will be applied, in case solutionMethod is "appliedMoment"
    balanceMomentNode = hingeNode+1

    # Flag for hinge rotation magnitude being known (fixed)
    rotationIsFixed = !isnothing(pHValue)

    # Normalize local hinge axis
    localHingeAxis /= norm(localHingeAxis)

    # Set initial (undeformed beam) and deformed hinge axes resolved in basis A
    deformedHingeAxis = initialHingeAxis = beam.R0*localHingeAxis

    # Set master (free) and slave directions in case of unknonwn hinge rotation
    slaveDir = [1,2,3]
    if rotationIsFixed
        masterDir = 0
    else
        masterDir = argmax(abs.(initialHingeAxis))
        popat!(slaveDir,masterDir)
    end

    # Initialize global IDs (updated later on model assembly)
    masterGlobalID, slaveGlobalID = 0, 0

    # Initialize global DOFs indices (updated later on model assembly)
    masterGlobalDOFs, slaveGlobalDOFs, constraintGlobalDOFs, λEqs = zeros(Int,3), zeros(Int,3), zeros(Int,6), zeros(Int,length(slaveDir))

    # Initialize the balance moment vector (moment necessary to enforce the constraint)
    balanceMoment = rotationIsFixed ? zeros(3) : zeros(2)
    balanceMomentDirections = rotationIsFixed ? ["M1A","M2A","M3A"] : ifelse(masterDir==1,["M2A","M3A"],ifelse(masterDir==2,["M1A","M3A"],["M1A","M2A"]))

    # Initialize ID of the BC corresponding to the balance moment
    balanceMomentBCID = 0

    # Initialize the array of Lagrange multipliers
    λ = balanceMoment

    # Initialize Jacobian matrix of the constraint (size is updated later on model assembly)
    Jc = zeros(0,0)

    # Initialize vector of rotation parameters across the hinge, resolved in basis A, angle of rotation [rad] and hinge moment
    pH, ϕ, hingeMoment = zeros(3), 0.0, 0.0

    return HingeAxisConstraint(beam=beam,masterOnLeft=masterOnLeft,localHingeAxis=localHingeAxis,pHValue=pHValue,solutionMethod=solutionMethod,updateAllDOFinResidual=updateAllDOFinResidual,rotationIsFixed=rotationIsFixed,initialHingeAxis=initialHingeAxis,deformedHingeAxis=deformedHingeAxis,masterLocalID=masterLocalID,slaveLocalID=slaveLocalID,masterDir=masterDir,slaveDir=slaveDir,masterGlobalID=masterGlobalID, slaveGlobalID=slaveGlobalID,masterGlobalDOFs=masterGlobalDOFs,slaveGlobalDOFs=slaveGlobalDOFs,constraintGlobalDOFs=constraintGlobalDOFs,λEqs=λEqs,balanceMomentNode=balanceMomentNode,balanceMoment=balanceMoment,balanceMomentDirections=balanceMomentDirections,balanceMomentBCID=balanceMomentBCID,λ=λ,Jc=Jc,pH=pH,ϕ=ϕ,hingeMoment=hingeMoment)

end
export create_HingeAxisConstraint


#
# mutable struct HingeAxisConstraintData

#     HingeAxisConstraintData composite type

# Fields
# `balanceMoment::Vector{Float64}`: balance moment vector around the hinge, resolved in basis A
# `hingeMoment::Real`: total moment about the hinge axis
# `pH::Vector{Float64}`: hinge rotation parameters, resolved in basis A
# `ϕ::Real`: hinge angle
#
mutable struct HingeAxisConstraintData
    
    # Fields
    balanceMoment::Vector{Float64}
    hingeMoment::Real
    pH::Vector{Float64}
    ϕ::Real

    # Constructor
    function HingeAxisConstraintData(balanceMoment,hingeMoment,pH,ϕ)
        return new(balanceMoment,hingeMoment,pH,ϕ)
    end
end