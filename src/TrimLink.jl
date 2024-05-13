"""
@with_kw mutable struct TrimLoadsLink

TrimLoadsLink composite type

# Fields
- masterBC::BC
- slaveBCs::Vector{BC}
- masterBeam::Beam
- masterNodeLocalID::Int64
- slaveBeams::Vector{Beam}
- slaveNodesLocalIDs::Vector{Int64}
- masterNodeGlobalID::Int64
- slaveNodesGlobalIDs::Vector{Int64}
"""
@with_kw mutable struct TrimLoadsLink

    # Primary fields (inputs)
    masterBC::BC
    slaveBCs::Vector{BC}
    # Secondary fields (outputs)
    masterBeam::Beam
    masterNodeLocalID::Int64
    slaveBeams::Vector{Beam}
    slaveNodesLocalIDs::Vector{Int64}
    masterNodeGlobalID::Int64
    slaveNodesGlobalIDs::Vector{Int64}

end

# Constructor
function create_TrimLoadsLink(;masterBC::Union{Nothing,BC},slaveBCs::Union{Nothing,Vector{BC}})

    # Get BCs information
    masterBeam = masterBC.beam
    masterNodeLocalID = masterBC.node
    slaveBeams = [BC.beam for BC in slaveBCs]
    slaveNodesLocalIDs = [BC.node for BC in slaveBCs]

    # Validate
    validLoadTypes = ["F1A","F2A","F3A","M1A","M2A","M3A","Ff1A","Ff2A","Ff3A","Mf1A","Mf2A","Mf3A"]
    for masterLoadType in masterBC.types
        @assert masterLoadType in validLoadTypes "invalid load type for master load BC"
    end
    for slaveBC in slaveBCs
        @assert length(slaveBC.types) == length(masterBC.types) "set one slave BC type for each master BC type"
        for slaveLoadType in slaveBC.types
            @assert slaveLoadType in validLoadTypes "invalid load type for slave load BC"
        end
    end
    @assert any(masterBC.isTrim) "master load BC does not contain any trim loads"
    for slaveBC in slaveBCs
        @assert any(slaveBC.isTrim) "slave load BC does not contain any trim loads"
    end

    # Global IDs are updated upon assembly of the model

    return TrimLoadsLink(masterBC,slaveBCs,masterBeam,masterNodeLocalID,slaveBeams,slaveNodesLocalIDs,0,Vector{Int64}())

end
export create_TrimLoadsLink