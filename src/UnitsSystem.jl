"""
@with_kw mutable struct UnitsSystem 

A composite type with fields for length, force, angle and frequency units (this is only for plotting purposes and does not influence calculations)

# Fields:
- length::String
- force::String
- angle::String
- frequency::String
"""
@with_kw mutable struct UnitsSystem

    # Fields
    length::String = "m"
    force::String = "N"
    angle::String = "rad"
    frequency::String = "rad/s"

    # Constructor
    function UnitsSystem(length::String,force::String,angle::String,frequency::String)

        self = new(length,force,angle,frequency)

        # Validate
        validate_unit_system(self)

        return self

    end
end


"""
validate_unit_system(units::UnitsSystem)

Validates the unit system

# Fields:
- units::UnitsSystem
"""
function validate_unit_system(units::UnitsSystem)

    @unpack length,force,angle,frequency = units

    # Length 
    possible_length_units = ["m","cm","mm","ft","in"]
    if !(length in possible_length_units)
        error("'length' must be one of $(possible_lengths)")
    end 
    
    # Force 
    possible_force_units = ["N","lbf"]
    if !(force in possible_force_units)
        error("'force' must be one of $(possible_force_units)")
    end

    # Angle 
    possible_angle_units = ["rad","deg"]
    if !(angle in possible_angle_units)
        error("'angle' must be one of $(possible_angle_units)")
    end

    # Frequency 
    possible_frequency_units = ["rad/s","Hz","rpm"]
    if !(frequency in possible_frequency_units)
        error("'frequency' must be one of $(possible_frequency_units)")
    end

end
