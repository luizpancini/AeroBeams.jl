## LINKS COMPOSITE TYPES
# ------------------------------------------------------------------------------
# TrimLoadsLink master and slave composite type
# Fields
# - `masterBC`: master BC
# - `masterBeam`: master beam
# - `slaveBCs`: slave BCs
# - `slaveBeams`: slave beams
#
struct TrimLoadsLinkMasterAndSlave

    masterBC::BC
    masterBeam::Beam

    slaveBCs::Vector{BC}
    slaveBeams::Vector{Beam}

end


# TrimLoadsLink assembly data composite type
# Fields
# - `masterNodeGlobalID`: global ID of the master node
# - `slaveNodesGlobalIDs`: global IDs of the slave nodes
#
struct TrimLoadsLinkAssemblyData

    masterNodeGlobalID::Int
    slaveNodesGlobalIDs::Vector{Int}

end
TrimLoadsLinkAssemblyData() = TrimLoadsLinkAssemblyData(
    0, zeros(Int,0)
)


# TrimLoadsLink composite type
# Fields
# - `masterAndSlave`: master and slave data for the trim loads link
# - `assemblyData`: assembly data for the trim loads link
#
@with_kw mutable struct TrimLoadsLink

    # Definition fields
    masterAndSlave::TrimLoadsLinkMasterAndSlave

    # Assembly-dependent field (to be set by Model)
    assemblyData::TrimLoadsLinkAssemblyData = TrimLoadsLinkAssemblyData()

end


# FlapLink composite type
# Fields
# - `masterBeam`: master beam
# - `slaveBeams`: slave beams
# - `δMultipliers`: multiplication factors of flap deflection in slave surfaces relative to the master surface
#
@with_kw struct FlapLink

    # Definition fields (inputs for FlapLink creation)
    masterBeam::Beam
    slaveBeams::Vector{Beam}
    δMultipliers::Vector{Float64}

end


## PUBLIC API
# ------------------------------------------------------------------------------
"""
    create_TrimLoadsLink(; kwargs...)

Trim loads link constructor (a link between trim loads so that they are equal)

# Keyword arguments
- `masterBC::BC`: master BC
- `slaveBCs::Vector{BC}`: slave BCs
"""
function create_TrimLoadsLink(; 
    masterBC::BC,
    slaveBCs::Vector{BC}
    )

    # Validate
    validLoadTypes = [:F1A,:F2A,:F3A,:M1A,:M2A,:M3A,:F1b,:F2b,:F3b,:M1b,:M2b,:M3b,:Ff1A,:Ff2A,:Ff3A,:Mf1A,:Mf2A,:Mf3A,:Ff1b,:Ff2b,:Ff3b,:Mf1b,:Mf2b,:Mf3b]
    for masterLoadType in masterBC.definition.types
        @assert masterLoadType in validLoadTypes "invalid load type for masterBC, must be one of: $(validLoadTypes)"
    end
    for (i, slaveBC) in enumerate(slaveBCs)
        for (j, type) in enumerate(slaveBC.definition.types)
            @assert type == masterBC.definition.types[j] "slaveBC $(i) types do not match masterBC's"
        end
    end
    @assert any(masterBC.deducedData.isTrim) "masterBC does not contain any trim loads"
    for (i, slaveBC) in enumerate(slaveBCs)
        @assert slaveBC.deducedData.isTrim == masterBC.deducedData.isTrim "slaveBC $(i) isTrim flags do not match masterBC's"
    end

    # Get beams from BCs
    masterBeam = masterBC.definition.beam
    slaveBeams = [BC.definition.beam for BC in slaveBCs]

    # Set master and slave data
    masterAndSlave = TrimLoadsLinkMasterAndSlave(masterBC,masterBeam,slaveBCs,slaveBeams)

    return TrimLoadsLink(
        masterAndSlave=masterAndSlave,
        )

end
export create_TrimLoadsLink


"""
    create_FlapLink(; kwargs...)

Flap link constructor (a link between flapped surfaces)

# Keyword arguments
- `masterBeam::Beam`: beam of the master surface
- `slaveBeams::Vector{Beam}`: beams of the slave surfaces
- `δMultipliers::Vector{<:Real}`: multiplication factors of flap deflection in slave surfaces relative to the master surface
"""
function create_FlapLink(; 
    masterBeam::Beam,
    slaveBeams::Vector{Beam},
    δMultipliers::Vector{<:Real}=ones(length(slaveBeams))
    )

    # Validate 
    @assert !isnothing(masterBeam.attachments.surface) "master beam has no aerodynamic surface"
    @assert !isnothing(masterBeam.attachments.surface.geometry.normFlapSpan) "master beam's aerodynamic surface has no flap"
    for beam in slaveBeams
        @assert beam !== masterBeam "beam set as both slave and master"
        @assert !isnothing(beam.attachments.surface) "slave beam has no aerodynamic surface"
        @assert !isnothing(beam.attachments.surface.geometry.normFlapSpan) "slave beam's aerodynamic surface has no flap"
    end
    @assert length(δMultipliers) == length(slaveBeams) "specify one deflection multiplier for each slave surface"

    # Set correct type
    δMultipliers = Float64.(δMultipliers)

    return FlapLink(
        masterBeam,
        slaveBeams,
        δMultipliers
        )

end
export create_FlapLink