#

#     RotationConstraint composite type

#
@with_kw mutable struct RotationConstraint

    # Primary (inputs to rotation constraint creation)
    beam::Beam
    masterElementLocalID::Int64
    slaveElementLocalID::Int64
    masterDOF::Int64
    slaveDOF::Int64
    value::Real
    loadBalanceLocalNode::Int64

    # Secondary (outputs from rotation constraint creation)
    masterElementGlobalID::Int64
    slaveElementGlobalID::Int64
    masterElemMasterGlobalDOF::Int64
    slaveElemSlaveGlobalDOF::Int64
    balanceMoment::Real

end
export RotationConstraint


"""
    create_RotationConstraint(; kwargs...)

Rotation constraint constructor

# Keyword arguments
- `beam::Beam` = beam
- `masterElementLocalID::Int64` = local ID of the master element
- `slaveElementLocalID::Int64` = local ID of the slave element
- `masterDOF::Int64` = master (reference) degree-of-freedom
- `slaveDOF::Int64` = slave degree-of-freedom
- `value::Real` = value of the slave DOF relative to the master DOF (Wiener-Milenkovic parameter value)
- `loadBalanceLocalNode::Int64` = local node of the beam where to set the balance load (i.e., the moment necessary to enforce the constraint)
"""
function create_RotationConstraint(; beam::Beam,masterElementLocalID::Int64,slaveElementLocalID::Int64,masterDOF::Int64,slaveDOF::Int64,value::Real,loadBalanceLocalNode::Int64)

    # Validate
    @assert masterElementLocalID <= beam.nElements "beam does not have 'masterElementLocalID' elements"
    @assert slaveElementLocalID <= beam.nElements "beam does not have 'slaveElementLocalID' elements"
    @assert 1 <= masterDOF <= 3 "specify masterDOF between 1 and 3"
    @assert 1 <= slaveDOF <= 3 "specify slaveDOF between 1 and 3"
    @assert -4 < value < 4 "set value between -4 and 4 (equivalent to rotation between -π and π)"

    # Initialize global IDs of elements (updated later on model assembly)
    masterElementGlobalID,slaveElementGlobalID = 0,0

    # Initialize global DOFs (updated later on model assembly)
    masterElemMasterGlobalDOF,slaveElemSlaveGlobalDOF = 0,0

    # Initialize the balance load (moment necessary to enforce the constraint)
    balanceMoment = 0.0

    return RotationConstraint(beam=beam,masterElementLocalID=masterElementLocalID,slaveElementLocalID=slaveElementLocalID,masterDOF=masterDOF,slaveDOF=slaveDOF,value=value,loadBalanceLocalNode=loadBalanceLocalNode,masterElementGlobalID=masterElementGlobalID,slaveElementGlobalID=slaveElementGlobalID,masterElemMasterGlobalDOF=masterElemMasterGlobalDOF,slaveElemSlaveGlobalDOF=slaveElemSlaveGlobalDOF,balanceMoment=balanceMoment)

end
export create_RotationConstraint