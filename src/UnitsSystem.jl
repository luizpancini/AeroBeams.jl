#
# @with_kw mutable struct UnitsSystem 

#     UnitsSystem composite type

#
@with_kw mutable struct UnitsSystem

    # Fields
    length::String
    force::String
    angle::String
    frequency::String
    mass::String

end


"""
create_UnitsSystem(; kwargs...)

Creates a system composed of length, force, angle, frequency and mass units (this is only for plotting purposes and does not influence calculations)

# Keyword arguments
- `length::String` = length unit
- `force::String` = force unit
- `angle::String` = angle unit
- `frequency::String` = frequency unit
- `mass::String` = mass unit
"""
function create_UnitsSystem(;length::String="m",force::String="N",angle::String="rad",frequency::String="rad/s",mass::String="kg")

    # Initialize
    self = UnitsSystem(length=length,force=force,angle=angle,frequency=frequency,mass=mass)

    # Validate
    validate_units_system(self)

    return self

end
export create_UnitsSystem


# Validates the units system
function validate_units_system(units::UnitsSystem)

    @unpack length,force,angle,frequency,mass = units

    # Validate length 
    availableLengthUnits = ["m","cm","mm","ft","in"]
    if !(length in availableLengthUnits)
        error("'length' must be one of $(available_lengths)")
    end 

    # Validate force 
    availableForceUnits = ["N","lbf","kN","kip"]
    if !(force in availableForceUnits)
        error("'force' must be one of $(availableForceUnits)")
    end

    # Validate angle 
    availableAngleUnits = ["rad","deg"]
    if !(angle in availableAngleUnits)
        error("'angle' must be one of $(availableAngleUnits)")
    end

    # Validate frequency 
    availableFrequencyUnits = ["rad/s","Hz","rpm"]
    if !(frequency in availableFrequencyUnits)
        error("'frequency' must be one of $(availableFrequencyUnits)")
    end

    # Validate mass 
    availableMassUnits = ["kg","g","lb"]
    if !(mass in availableMassUnits)
        error("'mass' must be one of $(availableMassUnits)")
    end

end