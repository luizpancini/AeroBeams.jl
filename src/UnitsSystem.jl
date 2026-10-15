## UNITS SYSTEM COMPOSITE TYPE
# ------------------------------------------------------------------------------
# UnitsSystem composite type
struct UnitsSystem

    length::String
    force::String
    angle::String
    frequency::String
    mass::String

end


## PUBLIC API
# ------------------------------------------------------------------------------
"""
    create_UnitsSystem(; kwargs...)

Creates a system composed of length, force, angle, frequency and mass units (this is only for plotting purposes and does not influence calculations)

# Keyword arguments
- `length::String`: length unit
- `force::String`: force unit
- `angle::String`: angle unit
- `frequency::String`: frequency unit
- `mass::String`: mass unit
"""
function create_UnitsSystem(; 
    length::String="m",
    force::String="N",
    angle::String="rad",
    frequency::String="rad/s",
    mass::String="kg"
    )

    # Validate length 
    availableLengthUnits = ["m","cm","mm","ft","in"]
    if !(length in availableLengthUnits)
        error("'length' must be one of $(availableLengthUnits)")
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

    return UnitsSystem(
        length,
        force,
        angle,
        frequency,
        mass
        )

end
export create_UnitsSystem


"""
    feet2meters(x)

Converts a length from feet to meters

# Arguments
- `x`: length in feet
"""
feet2meters(x) = x * 0.3048
function feet2meters!(x)
    x .*= 0.3048
end
export feet2meters, feet2meters!


"""
    pounds2kilograms(x)

Converts a mass from pounds to kilograms

# Arguments
- `x`: mass in pounds
"""
pounds2kilograms(x) = x * 0.453592
function pounds2kilograms!(x)
    x .*= 0.453592
end
export pounds2kilograms, pounds2kilograms!


"""
    rads2Hz(x)

Converts a frequency from radians per second to Hertz

# Arguments
- `x`: frequency in radians per second
"""
rads2Hz(x) = x / (2π)
function rads2Hz!(x)
    x ./= 2π
end
export rads2Hz, rads2Hz!


"""
    Hz2rads(x)

Converts a frequency from Hertz to radians per second

# Arguments
- `x`: frequency in Hertz
"""
Hz2rads(x) = x * (2π)
function Hz2rads!(x)
    x .*= 2π
end
export Hz2rads, Hz2rads!


"""
    rpm2rads(x)

Converts a frequency from revolutions per minute to radians per second

# Arguments
- `x`: frequency in revolutions per minute
"""
rpm2rads(x) = x * (2π/60)
function rpm2rads!(x)
    x .*= 2π/60
end
export rpm2rads, rpm2rads!


"""
    rads2rpm(x)

Converts a frequency from radians per second to revolutions per minute

# Arguments
- `x`: frequency in radians per second
"""
rads2rpm(x) = x * (60/(2π))
function rads2rpm!(x)
    x .*= 60/(2π)
end
export rads2rpm, rads2rpm!