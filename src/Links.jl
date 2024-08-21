#
# @with_kw mutable struct TrimLoadsLink

# TrimLoadsLink composite type

#
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


"""
    create_TrimLoadsLink(; masterBC::BC,slaveBCs::Vector{BC})

Creates a link between trim loads so that they are equal

# Keyword arguments
- `masterBC::BC` = master BC
- `slaveBCs::Vector{BC}` = slave BCs
"""
function create_TrimLoadsLink(; masterBC::BC,slaveBCs::Vector{BC})

    # Get BCs information
    masterBeam = masterBC.beam
    masterNodeLocalID = masterBC.node
    slaveBeams = [BC.beam for BC in slaveBCs]
    slaveNodesLocalIDs = [BC.node for BC in slaveBCs]

    # Validate
    validLoadTypes = ["F1A","F2A","F3A","M1A","M2A","M3A","F1b","F2b","F3b","M1b","M2b","M3b","Ff1A","Ff2A","Ff3A","Mf1A","Mf2A","Mf3A","Ff1b","Ff2b","Ff3b","Mf1b","Mf2b","Mf3b"]
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


#
# @with_kw mutable struct FlapLink

#     FlapLink composite type
    
#
@with_kw mutable struct FlapLink

    masterBeam::Beam
    slaveBeams::Vector{Beam}
    δMultipliers::Vector{<:Number}

end


"""
    create_FlapLink(; kwargs...)

Creates a link between flapped surfaces

# Keyword arguments
- `masterBeam::Beam` = beam of the master surface
- `slaveBeams::Vector{Beam}` = beams of the slave surfaces
- `δMultipliers::Vector{<:Number}` = multiplication factors of flap deflection in slave surfaces relative to the master surface
"""
function create_FlapLink(; masterBeam::Beam,slaveBeams::Vector{Beam},δMultipliers::Vector{<:Number}=ones(length(slaveBeams)))

    # Validate 
    @assert !isnothing(masterBeam.aeroSurface) "master beam has no aerodynamic surface"
    @assert !isnothing(masterBeam.aeroSurface.normFlapSpan) "master beam's aerodynamic surface has no flap"
    for beam in slaveBeams
        @assert beam !== masterBeam "beam set as both slave and master"
        @assert !isnothing(beam.aeroSurface) "slave beam has no aerodynamic surface"
        @assert !isnothing(beam.aeroSurface.normFlapSpan) "slave beam's aerodynamic surface has no flap"
    end
    @assert length(δMultipliers) == length(slaveBeams) "specify one deflection multiplier for each slave surface"

    return FlapLink(masterBeam,slaveBeams,δMultipliers)

end
export create_FlapLink