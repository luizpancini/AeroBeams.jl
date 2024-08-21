#

#     RotationConstraint composite type

#
@with_kw mutable struct RotationConstraint

    # Primary (inputs to rotation constraint creation)
    masterBeam::Beam
    slaveBeam::Beam
    masterElementLocalID::Int64
    slaveElementLocalID::Int64
    DOF::Int64
    value::Number

    # Secondary (outputs from hinge creation)
    masterElementGlobalID::Int64
    slaveElementGlobalID::Int64

end


"""
    create_RotationConstraint(; kwargs...)

Rotation constraint constructor

# Keyword arguments
- `masterBeam::Beam` = master beam
- `slaveBeam::Beam` = slave beam
- `masterElementLocalID::Int64` = local ID of the master element
- `slaveElementLocalID::Int64` = local ID of the slave element
- `DOF::Int64` = constrained rotation degree-of-freedom
- `value::Number` = value of the slave DOF relative to the master DOF
"""
function create_RotationConstraint(; masterBeam::Beam,slaveBeam::Beam,masterElementLocalID::Int64,slaveElementLocalID::Int64,DOF::Int64,value::Number)

    # Validate
    @assert masterElementLocalID <= masterBeam.nElements "masterBeam does not have 'masterElementLocalID' elements"
    @assert slaveElementLocalID <= slaveBeam.nElements "slaveBeam does not have 'slaveElementLocalID' elements"
    @assert 1 <= DOF <= 3 "specify DOF between 1 and 3"
    @assert -3π/4 <= value <= 3π/4 "set value between -3π/4 and 3π/4"

    # Initialize global IDs of master and slave elements (updated later on model assembly)
    masterElementGlobalID = 0
    slaveElementGlobalID = 0

    return RotationConstraint(masterBeam=masterBeam,slaveBeam=slaveBeam,masterElementLocalID=masterElementLocalID,slaveElementLocalID=slaveElementLocalID,DOF=DOF,value=value,masterElementGlobalID=masterElementGlobalID,slaveElementGlobalID=slaveElementGlobalID)

end
export create_RotationConstraint