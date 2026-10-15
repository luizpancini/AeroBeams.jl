## ATMOSPHERE COMPOSITE TYPE
# ------------------------------------------------------------------------------
# Atmosphere composite type
# Fields:
# - `altitude::Float64`: altitude
# - `ρ::Float64`: air density
# - `μ::Float64`: dynamic viscosity
# - `a::Float64`: speed of sound
# - `g::Float64`: acceleration of gravity
#
@with_kw struct Atmosphere

    altitude::Float64 = 0.0
    ρ::Float64 = 1.225
    μ::Float64 = 1.7894e-5
    a::Float64 = 340.3
    g::Float64 = 9.80665

end
export Atmosphere


## PUBLIC API
# ------------------------------------------------------------------------------
"""
    standard_atmosphere(altitude::Real)

Atmosphere constructor based on the International Standard Atmosphere (https://en.wikipedia.org/wiki/International_Standard_Atmosphere)

# Arguments
- `altitude::Real`: altitude
"""
function standard_atmosphere(altitude::Real)

    # Check altitude
    @assert altitude >= 0 "altitude must be non-negative"
    @assert altitude <= 47e3 "cannot calculate ISA for an altitude greater than 47 km"

    # Mean sea-level values
    #---------------------------------------------------------------------------
    # Pressure [Pa]
    p₀ = 101325
    # Temperature [K]
    T₀ = 288.15
    # Acceleration of gravity [m/s^2]
    g₀ = 9.80665
    
    # Other constants
    #---------------------------------------------------------------------------
    # Specific gas constant of air [J/kg.K]
    R = 287.05864074988347
    # Specific heat ratio
    γ = 1.4
    # Constants from Sutherland's law
    β,S = 1.458e-6, 110.4
    # Troposphere lapse rate [K/m] 
    λTropo = -6.5e-3
    # Stratosphere lapse rate [K/m] 
    λStrat = 1e-3    
    # Stratopause lapse rate [K/m]
    λStratPause = 2.8e-3
    
    # Pressure and temperature at input altitude
    #---------------------------------------------------------------------------
    # Gravity [m/s^2] (Earth's radius ≈ 6.371e6 m)
    g = g₀ * (6.371e6/(6.371e6+altitude))^2
    # Temperature [K] and pressure [Pa] (the calculations already consider geopotential altitude - see https://en.wikipedia.org/wiki/International_Standard_Atmosphere)
    T₁₁ = T₀+λTropo*11e3
    T₂₀ = T₀-71.5
    T₃₂ = T₀-59.5
    p₁₁ = p₀*(1+λTropo*11e3/T₀)^(-g₀/(R*λTropo))
    p₂₀ = p₁₁*exp(-g₀/(R*T₁₁)*(20e3-11e3))
    p₃₂ = p₂₀*(1+λStrat*(32e3-20e3)/T₂₀)^(-g₀/(R*λStrat))
    if altitude <= 11e3
        T = T₀+λTropo*altitude
        p = p₀*(1+λTropo*altitude/T₀)^(-g₀/(R*λTropo))
    elseif altitude <= 20e3
        T = T₀-71.5 
        p = p₁₁*exp(-g₀/(R*T₁₁)*(altitude-11e3))
    elseif altitude <= 32e3
        T = T₀-71.5+λStrat*(altitude-20e3)
        p = p₂₀*(1+λStrat*(altitude-20e3)/T₂₀)^(-g₀/(R*λStrat))
    elseif altitude <= 47e3
        T = T₀-59.5+λStratPause*(altitude-32e3)
        p = p₃₂*(1+λStratPause*(altitude-32e3)/T₃₂)^(-g₀/(R*λStratPause))
    end
    
    # Corresponding air density, sound speed and dynamic viscosity
    ρ,μ,a = air_properties_from_pressure_and_temperature(p,T,R,γ,β,S)

    return Atmosphere(altitude=altitude,ρ=ρ,μ=μ,a=a,g=g)

end
export standard_atmosphere


"""
    air_properties_from_pressure_and_temperature(p::Real,T::Real,R::Real,γ::Real,β::Real,S::Real)

Calculates air properties from given pressure and temperature.

# Arguments
- `p::Real`: pressure
- `T::Real`: temperature
- `R::Real`: specific gas constant of air (default: 287.05864074988347)
- `γ::Real`: specific heat ratio (default: 1.4)
- `β::Real`: Sutherland's law constant (default: 1.458e-6)
- `S::Real`: Sutherland's law constant (default: 110.4)

# Returns
- `ρ::Float64`: air density
- `μ::Float64`: dynamic viscosity
- `a::Float64`: speed of sound
"""
function air_properties_from_pressure_and_temperature(p::Real,T::Real,R::Real=287.05864074988347,γ::Real=1.4,β::Real=1.458e-6,S::Real=110.4)

    # Check inputs
    @assert all(x -> x>0, (p,T,R,γ,β,S)) "all inputs must be positive real numbers"

    # Density
    ρ = p/(R*T)

    # Dynamic viscosity
    μ = β/(T+S)*T^(3/2)

    # Sound speed
    a = sqrt(γ*R*T)

    return ρ,μ,a

end
export air_properties_from_pressure_and_temperature