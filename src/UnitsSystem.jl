"""
@with_kw mutable struct UnitsSystem 

A composite type with fields for length, force, angle, frequency and mass units (this is only for plotting purposes and does not influence calculations)

# Fields:
- length::String
- force::String
- angle::String
- frequency::String
- mass::String
"""
@with_kw mutable struct UnitsSystem

    # Fields
    length::String
    force::String
    angle::String
    frequency::String
    mass::String

end

# Constructor
function create_UnitsSystem(;length::String="m",force::String="N",angle::String="rad",frequency::String="rad/s",mass::String="kg")

    # Initialize
    self = UnitsSystem(length=length,force=force,angle=angle,frequency=frequency,mass=mass)

    # Validate
    validate_units_system(self)

    return self

end
export create_UnitsSystem


"""
validate_units_system(units::UnitsSystem)

Validates the units system

# Arguments
- `units::UnitsSystem`
"""
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