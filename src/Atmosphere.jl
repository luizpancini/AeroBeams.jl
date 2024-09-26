#

    # Atmosphere composite type

#
@with_kw mutable struct Atmosphere

    ρ::Real = 1.225
    μ::Real = 1.7894e-5
    a::Real = 340.3

end
export Atmosphere


# Calculates air properties from given pressure and temperature
function air_properties_from_pressure_and_temperature(p::Real,T::Real,R::Real=287.05864074988347,γ::Real=1.4,β::Real=1.458e-6,S::Real=110.4)

    # Check inputs
    @assert p > 0
    @assert T > 0
    @assert R > 0
    @assert γ > 0
    @assert β > 0
    @assert S > 0

    # Density [kg/m^3]
    ρ = p/(R*T)

    # Dynamic viscosity [Pa.s] 
    μ = β/(T+S)*T^(3/2)

    # Sound speed [m/s]
    a = sqrt(γ*R*T)

    return ρ,μ,a

end


"""
    standard_atmosphere(altitude::Real)

Atmosphere constructor based on the International Standard Atmosphere (https://en.wikipedia.org/wiki/International_Standard_Atmosphere)

# Arguments
- `altitude::Real`
"""
function standard_atmosphere(altitude::Real)

    # Check altitude
    if altitude > 47e3
        error("Cannot calculate ISA for an altitude greater than 47 km")
    end

    # Mean-sea-level values
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
    # Gravity [m/s^2] (the calculations already consider geopotential altitude - see https://en.wikipedia.org/wiki/International_Standard_Atmosphere)
    g = g₀
    # Temperature [K] and pressure [Pa]
    T₁₁ = T₀+λTropo*11e3
    T₂₀ = T₀-71.5
    T₃₂ = T₀-59.5
    p₁₁ = p₀*(1+λTropo*11e3/T₀)^(-g/(R*λTropo))
    p₂₀ = p₁₁*exp(-g/(R*T₁₁)*(20e3-11e3))
    p₃₂ = p₂₀*(1+λStrat*(32e3-20e3)/T₂₀)^(-g/(R*λStrat))
    if altitude <= 11e3
        T = T₀+λTropo*altitude
        p = p₀*(1+λTropo*altitude/T₀)^(-g/(R*λTropo))
    elseif altitude <= 20e3
        T = T₀-71.5 
        p = p₁₁*exp(-g/(R*T₁₁)*(altitude-11e3))
    elseif altitude <= 32e3
        T = T₀-71.5+λStrat*(altitude-20e3)
        p = p₂₀*(1+λStrat*(altitude-20e3)/T₂₀)^(-g/(R*λStrat))
    elseif altitude <= 47e3
        T = T₀-59.5+λStratPause*(altitude-32e3)
        p = p₃₂*(1+λStratPause*(altitude-32e3)/T₃₂)^(-g/(R*λStratPause))
    end
    
    # Corresponding air density, sound speed and dynamic viscosity
    ρ,μ,a = air_properties_from_pressure_and_temperature(p,T)

    # Set atmosphere
    atmosphere = Atmosphere(ρ,μ,a)

    return atmosphere

end
export standard_atmosphere