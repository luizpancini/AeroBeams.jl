abstract type Gust end

#
# @with_kw mutable struct SharpEdgedGust <: Gust

#     SharpEdgedGust composite type

#
@with_kw mutable struct SharpEdgedGust <: Gust

    # Primary (necessary for gust creation)
    initialTime::Real
    duration::Real
    convectiveVelocity::Real
    verticalVelocity::Real
    p::Vector{<:Real}
    
    # Secondary (fixed or outputs of gust creation)
    isDefinedOverTime::Bool = true
    UGustInertial::Function
    finalTime::Real

end


"""
    create_SharpEdgedGust(; kwargs...)

Creates a sharp-edged gust

# Keyword arguments
- `initialTime::Real`: time when the gust begins
- `duration::Real`: duration of the gust
- `convectiveVelocity::Real`: convective velocity of the gust (in longitudinal direction)
- `verticalVelocity::Real`: peak vertical velocity of the gust (in lift direction)
- `p::Vector{<:Real}`: Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
"""
function create_SharpEdgedGust(; initialTime::Real=0,duration::Real=Inf,convectiveVelocity::Real=0,verticalVelocity::Real,p::Vector{<:Real}=zeros(3))

    # Validate
    @assert duration > 0
    @assert initialTime >= 0
    @assert size(p) == (3,)

    # Rotation tensor from the inertial basis to the gust basis
    R = rotation_tensor_E321(p)
    RT = R'

    # Final time of the gust
    finalTime = initialTime + duration
    
    # Set vector of gust velocity (in the inertial frame) function over time
    UGustInertial = t::Real -> initialTime<t<finalTime ? RT*[0; convectiveVelocity; verticalVelocity] : zeros(3)

    return SharpEdgedGust(initialTime=initialTime,duration=duration,convectiveVelocity=convectiveVelocity,verticalVelocity=verticalVelocity,p=p,finalTime=finalTime,UGustInertial=UGustInertial)

end
export create_SharpEdgedGust


#
# @with_kw mutable struct OneMinusCosineGust <: Gust

#     OneMinusCosineGust composite type

#
@with_kw mutable struct OneMinusCosineGust <: Gust

    # Primary (necessary for gust creation)
    initialTime::Real
    duration::Real
    convectiveVelocity::Real
    verticalVelocity::Real
    p::Vector{<:Real}
    
    # Secondary (fixed or outputs of gust creation)
    isDefinedOverTime::Bool = true
    UGustInertial::Function
    finalTime::Real

end


"""
    create_OneMinusCosineGust(; kwargs...)

Creates a one-minus-cosine gust

# Keyword arguments
- `initialTime::Real`: time when the gust begins
- `duration::Real`: duration of the gust
- `convectiveVelocity::Real`: convective velocity of the gust (in longitudinal direction)
- `verticalVelocity::Real`: peak vertical velocity of the gust (in lift direction)
- `p::Vector{<:Real}`: Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
"""
function create_OneMinusCosineGust(; initialTime::Real=0,duration::Real,convectiveVelocity::Real=0,verticalVelocity::Real,p::Vector{<:Real}=zeros(3))

    # Validate
    @assert duration > 0
    @assert initialTime >= 0
    @assert size(p) == (3,)

    # Rotation tensor from the inertial basis to the gust basis
    R = rotation_tensor_E321(p)
    RT = R'

    # Final time of the gust
    finalTime = initialTime + duration
    
    # Set vector of gust velocity (in the inertial frame) function over time
    UGustInertial = t::Real -> initialTime<t<finalTime ? RT*[0; convectiveVelocity; 1/2*verticalVelocity*(1-cos(2π*(t-initialTime)/duration))] : zeros(3)

    return OneMinusCosineGust(initialTime=initialTime,duration=duration,convectiveVelocity=convectiveVelocity,verticalVelocity=verticalVelocity,p=p,finalTime=finalTime,UGustInertial=UGustInertial)

end
export create_OneMinusCosineGust


#
# @with_kw mutable struct Continuous1DGust <: Gust

#     Continuous1DGust composite type

#
@with_kw mutable struct Continuous1DGust <: Gust

    # Primary (necessary for gust creation)
    spectrum::String
    generationMethod::String
    initialTime::Real
    duration::Real
    generationDuration::Real
    components::Vector{Int}
    L::Real
    ωmin::Real
    ωmax::Real
    Uref::Real
    convectiveVelocity::Real
    σ::Real
    p::Vector{<:Real}
    seed::Int
    plotPSD::Bool
    
    # Secondary (fixed or outputs of gust creation)
    isDefinedOverTime::Bool = true
    UGustInertial::Function
    finalTime::Real
    U::Function
    V::Function
    W::Function

end


"""
    create_Continuous1DGust(; kwargs...)

Creates a continuous 1D gust

# Keyword arguments
- `spectrum::String`: spectrum of the PSD: von Kármán ("vK") or Dryden ("Dryden")
- `generationMethod::String`: method for generation of gust velocity ("sinusoids" or "whiteNoise")
- `initialTime::Real`: time when the gust begins
- `duration::Real`: duration of the gust
- `generationDuration::Real`: duration of the gust considered for its generation
- `components::Vector{Int}`: components of gust velocity to consider
- `L::Real`: turbulence length scale [m]
- `ωmin::Real`: minimum frequency of the PSD [rad/s]
- `ωmax::Real`: maximum frequency of the PSD [rad/s]
- `Uref::Real`: reference airspeed
- `convectiveVelocity::Real`: convective velocity of the gust (in longitudinal direction)
- `σ::Real`: turbulence intensity (RMS) of gust velocity component
- `p::Vector{<:Real}`: Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
- `seed::Int`: seed for random numbers generation (reproducibility)
- `plotPSD::Bool`: flag to plot the PSD
"""
function create_Continuous1DGust(; spectrum::String="vK",generationMethod::String="sinusoids",initialTime::Real=0,duration::Real,generationDuration::Real=duration,components::Vector{Int}=[3],L::Real=762,ωmin::Real=0,ωmax::Real=200π,Uref::Real,convectiveVelocity::Real=0,σ::Real,p::Vector{<:Real}=zeros(3),seed::Int=123456,plotPSD::Bool=false)

    # Validate
    @assert spectrum in ["vK","Dryden"]
    @assert generationMethod in ["sinusoids","whiteNoise"]
    @assert components in [[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]
    @assert initialTime >= 0
    @assert 0 < duration <= generationDuration
    @assert L > 0
    @assert 0 <= ωmin < ωmax
    @assert Uref > 0
    @assert σ > 0
    @assert size(p) == (3,)
    @assert seed > 0

    # Print warning
    if generationDuration > duration
        println("Gust duration implied for its generation is greater than actual duration: gust intensity and PSD might not correspond to desired input")
    end

    # Set random seed
    Random.seed!(seed)

    # Get vertical velocity according to generation method
    if generationMethod == "sinusoids"
        # Frequency resolution [rad/s]: this ensures a stable value on the RMS of the gust velocity profile
        dω = min(2π/generationDuration, (2π)^2/ωmax) 
        # Frequency vector [rad/s]
        ω = collect(ωmin:dω:ωmax)
        # Spatial frequency vector and resolution [rad/m]
        Ω = ω/Uref
        dΩ = dω/Uref
        # PSD [(m/s)^2/(rad/m)]
        if spectrum == "vK"
            Φu = @. σ^2*2*L/π * 1 / ( 1+(1.339*L*Ω)^2 )^(5/6)
            Φv = Φw = @. σ^2*L/π * ( (1+8/3*(1.339*L*Ω)^2) / (1+(1.339*L*Ω)^2)^(11/6) )
        elseif spectrum == "Dryden"
            Φu = @. σ^2*2*L/π * 1 / (1+(L*Ω)^2)^2
            Φv = Φw = @. σ^2*L/π * ( (1+3*(L*Ω)^2) / (1+(L*Ω)^2)^2 )
        end
        sqrt2dΩΦu = sqrt.(2*dΩ*Φu)
        sqrt2dΩΦv = sqrt.(2*dΩ*Φv)
        sqrt2dΩΦw = sqrt.(2*dΩ*Φw)
        # Random phase vectors
        ψu = 2π*rand(length(ω))
        ψv = 2π*rand(length(ω))
        ψw = 2π*rand(length(ω))
        # Gust velocity components as functions of time
        U = 1 in components ? t -> sum( sqrt2dΩΦu .* cos.(ω*t.+ψu) ) : t -> 0
        V = 2 in components ? t -> sum( sqrt2dΩΦv .* cos.(ω*t.+ψv) ) : t -> 0
        W = 3 in components ? t -> sum( sqrt2dΩΦw .* cos.(ω*t.+ψw) ) : t -> 0
        # Plot PSD of normal component, if applicable
        if plotPSD
            Δt = π/ωmax
            time = collect(0:Δt:generationDuration-Δt)
            f,_,UPSD = get_FFT_and_PSD(time,U.(time))
            f,_,WPSD = get_FFT_and_PSD(time,W.(time))
            plt = plot(xscale=:log10,yscale=:log10,xlims=[L*Ω[2],L*Ω[end]],xticks=[1e0,1e1,1e2,1e3,1e4,1e5,1e6],xlabel="Normalized spatial frequency (\$\\Omega L\$)",ylabel="Normalized PSD",tickfont=font(10),guidefont=font(16),legendfontsize=12)
            scatter!(L*Ω,Φu/(σ^2*2*L/π),c=:green,ms=5,msw=0,label="Definition - \$\\Phi_u\$")
            scatter!(L*Ω,Φw/(σ^2*L/π),c=:blue,ms=5,msw=0,label="Definition - \$\\Phi_w\$")
            plot!(L*Ω,UPSD/(2π)/(σ^2*2*L/(π*Uref)),c=:green,lw=2,label="Generated - \$\\Phi_u\$")
            plot!(L*Ω,WPSD/(2π)/(σ^2*L/(π*Uref)),c=:blue,lw=2,label="Generated - \$\\Phi_w\$")
            display(plt)
        end
    elseif generationMethod == "whiteNoise"
        # Time vector 
        Δt = π/ωmax
        time = collect(0:Δt:generationDuration-Δt)
        # Gust velocity arrays and PSDs
        Uwn,Vwn,Wwn,f,Φu,Φv,Φw = stochastic_gust_velocity_from_white_noise(spectrum,time,Uref,L,σ)
        # Mapping functions
        itpU = Interpolate(time, Uwn)
        itpV = Interpolate(time, Vwn)
        itpW = Interpolate(time, Wwn)
        # Gust velocity components as functions of time
        U = 1 in components ? t -> itpU(t) : t -> 0
        V = 2 in components ? t -> itpV(t) .+ convectiveVelocity : t -> 0
        W = 3 in components ? t -> itpW(t) : t -> 0
        # Plot PSD of normal gust velocity, if applicable
        if plotPSD
            f2,_,UPSD = get_FFT_and_PSD(time,Uwn)
            f2,_,WPSD = get_FFT_and_PSD(time,Wwn)
            LΩ = 2π*f/Uref*L
            LΩ2 = 2π*f2/Uref*L
            plt = plot(xscale=:log10,yscale=:log10,xlims=[LΩ[2],LΩ[end]],xlabel="Normalized spatial frequency",ylabel="Normalized PSD",tickfont=font(10),guidefont=font(16),legendfontsize=12,legend=:bottomleft)
            plot!(LΩ2,UPSD/(σ^2*L/(π*Uref)),lw=2,color=:green,label="Generated - \$\\Phi_u\$")
            plot!(LΩ2,WPSD/(σ^2*L/(π*Uref)),lw=2,color=:blue,label="Generated - \$\\Phi_w\$")
            scatter!(LΩ,Φu/(σ^2*L/(π*Uref)),lw=2,ms=5,msw=0,color=:green,label="Definition - \$\\Phi_u\$")
            scatter!(LΩ,Φw/(σ^2*L/(π*Uref)),lw=2,ms=5,msw=0,color=:blue,label="Definition - \$\\Phi_w\$")
            display(plt)
        end
    end

    # Rotation tensor from the inertial basis to the gust basis
    R = rotation_tensor_E321(p)
    RT = R'

    # Final time of the gust
    finalTime = initialTime + duration

    # Set vector of gust velocity (in the inertial frame) function over time
    UGustInertial = t::Real -> initialTime<t<finalTime ? RT*[U(t); V(t); W(t)] : zeros(3)

    return Continuous1DGust(spectrum=spectrum,generationMethod=generationMethod,initialTime=initialTime,duration=duration,generationDuration=generationDuration,components=components,L=L,ωmin=ωmin,ωmax=ωmax,Uref=Uref,convectiveVelocity=convectiveVelocity,σ=σ,p=p,seed=seed,plotPSD=plotPSD,finalTime=finalTime,UGustInertial=UGustInertial,U=U,V=V,W=W)

end
export create_Continuous1DGust


#
# @with_kw mutable struct DiscreteSpaceGust <: Gust

#     DiscreteSpaceGust composite type

#
@with_kw mutable struct DiscreteSpaceGust <: Gust

    # Primary (necessary for gust creation)
    type::String
    gustLength::Real
    gustWidth::Real
    convectiveVelocity::Real
    verticalVelocity::Real
    c0::Vector{<:Real}
    p::Vector{<:Real}
    
    # Secondary (fixed or outputs of gust creation)
    isDefinedOverTime::Bool = false
    UGustInertial::Function

end


"""
    create_DiscreteSpaceGust(; kwargs...)

Creates a discrete space gust

# Keyword arguments
- `type::String`: type of gust
- `gustLength::Real`: length of the gust (in longitudinal direction)
- `gustWidth::Real`: width of the gust (in spanwise direction)
- `convectiveVelocity::Real`: convective velocity of the gust (in longitudinal direction)
- `verticalVelocity::Real`: peak vertical velocity of the gust (in lift direction)
- `c0::Vector{<:Real}`: position vector of the front of the gust, resolved in the inertial basis
- `p::Vector{<:Real}`: Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
"""
function create_DiscreteSpaceGust(; type::String,gustLength::Real,gustWidth::Real,convectiveVelocity::Real=0,verticalVelocity::Real,c0::Vector{<:Real}=zeros(3),p::Vector{<:Real}=zeros(3))

    # Validate
    @assert type in ["SharpEdged","OneMinusCosine","DARPA"]
    @assert gustLength > 0
    @assert gustWidth > 0
    @assert size(c0) == (3,)
    @assert size(p) == (3,)

    # Rotation tensor from the inertial basis to the gust basis
    R = rotation_tensor_E321(p)
    RT = R'

    # Normalized width coordinate (ranges from -1/2 to 1/2 for -width/2 <= x[1] <= width/2)
    w = x::AbstractVector{<:Real} -> (x[1]-c0[1])/gustWidth

    # Normalized length coordinate (ranges from 0 to 1 for 0 <= x[2] <= length)
    l = x::AbstractVector{<:Real} -> (x[2]-c0[2])/gustLength

    # Gust front longitudinal coordinate as a function of time
    c(t::Real) = c0[2] + convectiveVelocity*t

    # TF function for a position vector being inside the gust at a given time
    isInside(x::AbstractVector{<:Real},t::Real) = -gustWidth/2 <= (x[1]-c0[1]) <= gustWidth/2 && 0 <= x[2]-c(t) <= gustLength

    # Vector transformation from inertial to gust basis
    I2g = x::AbstractVector{<:Real} -> RT * x
    
    # Set vertical velocity function according to gust type
    if type == "SharpEdged"
        W = x::AbstractVector{<:Real} -> verticalVelocity
    elseif type == "OneMinusCosine"
        W = x::AbstractVector{<:Real} -> 1/2*verticalVelocity*(1-cos(2π*l(I2g(x))))
    elseif type == "DARPA"
        W = x::AbstractVector{<:Real} -> -1/2*verticalVelocity*cos(2π*w(I2g(x)))*(1-cos(2π*l(I2g(x))))
    end

    # Set vector of gust velocity (in the inertial frame) as a function of position and time
    UGustInertial(x::AbstractVector{<:Real},t::Real) = isInside(I2g(x),t) ? RT*[0; convectiveVelocity; W(x)] : zeros(3)

    return DiscreteSpaceGust(type=type,gustLength=gustLength,gustWidth=gustWidth,convectiveVelocity=convectiveVelocity,verticalVelocity=verticalVelocity,c0=c0,p=p,UGustInertial=UGustInertial)

end
export create_DiscreteSpaceGust


#
# @with_kw mutable struct Continuous1DSpaceGust <: Gust

#     Continuous1DSpaceGust composite type

#
@with_kw mutable struct Continuous1DSpaceGust <: Gust

    # Primary (necessary for gust creation)
    spectrum::String
    gustLength::Real
    N::Int
    L::Real
    σ::Real
    c0::Vector{<:Real}
    p::Vector{<:Real}
    seed::Int
    
    # Secondary (fixed or outputs of gust creation)
    isDefinedOverTime::Bool = false
    UGustInertial::Function
    U::Function
    V::Function
    W::Function

end


"""
    create_Continuous1DSpaceGust(; kwargs...)

Creates a continuous 1D space gust

# Keyword arguments
- `spectrum::String`: spectrum of the PSD: von Kármán ("vK") or Dryden ("Dryden")
- `gustLength::Real`: length of the gust (in longitudinal direction)
- `N::Int`: number of nodes for length discretization
- `L::Real`: turbulence length scale
- `σ::Real`: turbulence intensity (RMS) of gust velocity components
- `c0::Vector{<:Real}`: position vector of the front of the gust, resolved in the inertial basis
- `p::Vector{<:Real}`: Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
- `seed::Int`: seed for random numbers generation (reproducibility)
- `plotPSD::Bool`: flag to plot the PSD
"""
function create_Continuous1DSpaceGust(; spectrum::String,gustLength::Real,N::Int=1001,L::Real=762,σ::Real=1,c0::Vector{<:Real}=zeros(3),p::Vector{<:Real}=zeros(3),seed::Int=123456,plotPSD::Bool=false)

    # Validate
    @assert spectrum in ["vK","Dryden"]
    @assert gustLength > 0
    @assert N > 0
    @assert L > 0
    @assert size(c0) == (3,)
    @assert size(p) == (3,)

    # Rotation tensor from the inertial basis to the gust basis
    R = rotation_tensor_E321(p)
    RT = R'

    # TF function for a position vector being inside the gust
    isInside(x::AbstractVector{<:Real}) = 0 <= x[2]-c0[2] <= gustLength

    # Vector transformation from inertial to gust basis
    I2g = x::AbstractVector{<:Real} -> RT * x

    # Second (longitudinal) component of the above
    I2g2 = x::AbstractVector{<:Real} -> I2g(x)[2]

    # Non-dimensional spatial frequency vector
    LΩ = 2π*L/gustLength * LinRange(0,div(N-1,2),N)

    # 1D gust PSDs
    if spectrum == "vK"
        Φu = @. σ^2*2*L/π * 1 / ( 1+(1.339*LΩ)^2 )^(5/6)
        Φv = @. σ^2*L/π * ( (1+8/3*(1.339*LΩ)^2) / (1+(1.339*LΩ)^2)^(11/6) )
        Φw = Φv
    elseif spectrum == "Dryden"
        Φu = @. σ^2*2*L/π * 1 / (1+(LΩ)^2)^2
        Φv = @. σ^2*L/π * ( (1+3*(LΩ)^2) / (1+(LΩ)^2)^2 )
        Φw = Φv
    end

    # Set random seed
    Random.seed!(seed)

    # White noise signals
    un = randn(N)
    vn = randn(N)
    wn = randn(N)

    # Gust velocity components in the spatial frequency domain
    Un = fft(un) .* sqrt.(Φu)
    Vn = fft(vn) .* sqrt.(Φv)
    Wn = fft(wn) .* sqrt.(Φw)

    # Gust velocity components in the spatial domain (values at the nodes)
    ug = real(ifft(Un))
    vg = real(ifft(Vn))
    wg = real(ifft(Wn))

    # Re-scale
    ug = ug*σ/rms(ug)
    vg = vg*σ/rms(vg)
    wg = wg*σ/rms(wg)

    # Spatial coordinate 
    X = c0[2] .+ LinRange(0,gustLength,N)

    # Mapping functions
    itpU = Interpolate(X, ug)
    itpV = Interpolate(X, vg)
    itpW = Interpolate(X, wg)

    # Gust velocity components as functions of spatial coordinate
    U = x::AbstractVector{<:Real} -> itpU(I2g2(x))
    V = x::AbstractVector{<:Real} -> itpV(I2g2(x))
    W = x::AbstractVector{<:Real} -> itpW(I2g2(x))

    # Plot PSD, if applicable
    if plotPSD
        time = collect(LinRange(0,gustLength/(2π),N))
        Ω2,_,yPSD = get_FFT_and_PSD(time,wg,tol=1e-9)
        yPSDSmoothed = moving_average(yPSD,10)
        plt = plot(xscale=:log10,yscale=:log10,xlims=[LΩ[2],Inf],xticks=[1e0,1e1,1e2,1e3,1e4,1e5,1e6],xlabel="Normalized spatial frequency (\$\\Omega L\$)",ylabel="Normalized PSD",tickfont=font(10),guidefont=font(16),legendfontsize=12,legend=:bottomleft)
        plot!(Ω2*L,yPSD/(σ^2*L/π),lw=2,label="Generated")
        plot!(Ω2*L,yPSDSmoothed/(σ^2*L/π),lw=2,label="Generated - smoothed")
        plot!(LΩ,Φw/(σ^2*L/π),lw=2,label="Definition")
        display(plt)
    end

    # Set vector of gust velocity (in the inertial frame) as a function of position
    UGustInertial(x::AbstractVector{<:Real},t::Real) = isInside(I2g(x)) ? RT*[U(x); V(x); W(x)] : zeros(3)

    return Continuous1DSpaceGust(spectrum=spectrum,gustLength=gustLength,N=N,L=L,σ=σ,c0=c0,p=p,seed=seed,UGustInertial=UGustInertial,U=U,V=V,W=W)

end
export create_Continuous1DSpaceGust


#
# @with_kw mutable struct Continuous2DSpaceGust <: Gust

#     Continuous2DSpaceGust composite type

#
@with_kw mutable struct Continuous2DSpaceGust <: Gust

    # Primary (necessary for gust creation)
    spectrum::String
    gustLength::Real
    gustWidth::Real
    Nx::Int
    Ny::Int
    L::Real
    σ::Real
    c0::Vector{<:Real}
    p::Vector{<:Real}
    seed::Int
    
    # Secondary (fixed or outputs of gust creation)
    isDefinedOverTime::Bool = false
    UGustInertial::Function
    U::Function
    V::Function
    W::Function

end


"""
    create_Continuous2DSpaceGust(; kwargs...)

Creates a continuous 2D space gust

# Keyword arguments
- `spectrum::String`: spectrum of the PSD: von Kármán ("vK") or Dryden ("Dryden")
- `gustLength::Real`: length of the gust (in longitudinal direction)
- `gustWidth::Real`: width of the gust (in lateral direction)
- `Nx::Int`: number of nodes for length discretization
- `Ny::Int`: number of nodes for width discretization
- `L::Real`: turbulence length scale
- `σ::Real`: turbulence intensity (RMS) of gust velocity component
- `c0::Vector{<:Real}`: position vector of the front of the gust, resolved in the inertial basis
- `p::Vector{<:Real}`: Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
- `seed::Int`: seed for random numbers generation (reproducibility)
"""
function create_Continuous2DSpaceGust(; spectrum::String,gustLength::Real,gustWidth::Real,Nx::Int=101,Ny::Int=101,L::Real=762,σ::Real=1,c0::Vector{<:Real}=zeros(3),p::Vector{<:Real}=zeros(3),seed::Int=123456)

    # Validate
    @assert spectrum in ["vK","Dryden"]
    @assert gustLength > 0
    @assert gustWidth > 0
    @assert Nx > 0
    @assert Ny > 0
    @assert L > 0
    @assert size(c0) == (3,)
    @assert size(p) == (3,)

    # Rotation tensor from the inertial basis to the gust basis
    R = rotation_tensor_E321(p)
    RT = R'

    # TF function for a position vector being inside the gust
    isInside(x::AbstractVector{<:Real}) = -gustWidth/2 <= x[1]-c0[1] <= gustWidth/2 && 0 <= x[2]-c0[2] <= gustLength

    # Vector transformation from inertial to gust basis
    I2g = x::AbstractVector{<:Real} -> RT * x

    # Longitudinal and lateral components of the above
    I2g12 = x::AbstractVector{<:Real} -> I2g(x)[1:2]

    # Constant
    if spectrum == "vK"
        a = 1.339
    elseif spectrum == "Dryden"
        a = 1
    end

    # Non-dimensional spatial frequency vectors
    aLΩx = 2π*a*L/gustWidth * LinRange(0,div(Nx-1,2),Nx)
    aLΩy = 2π*a*L/gustLength * LinRange(0,div(Ny-1,2),Ny)

    # 2D gust PSD matrices
    if spectrum == "vK"
        # Uncorrelated PSDs (Etkin - Dynamics of Atmospheric Flight)
        Φuu = [(σ*a*L)^2/(6π)/(1+aLΩx[i]^2+aLΩy[j]^2)^(7/3) * (1+aLΩx[i]^2+11/3*aLΩy[j]^2) for i in eachindex(aLΩx), j in eachindex(aLΩy)]
        Φvv = [(σ*a*L)^2/(6π)/(1+aLΩx[i]^2+aLΩy[j]^2)^(7/3) * (1+11/3*aLΩx[i]^2+aLΩy[j]^2) for i in eachindex(aLΩx), j in eachindex(aLΩy)]
        Φww = [(σ*a*L)^2/(9π)/(1+aLΩx[i]^2+aLΩy[j]^2)^(7/3) * 4*(aLΩx[i]^2+aLΩy[j]^2) for i in eachindex(aLΩx), j in eachindex(aLΩy)]
    elseif spectrum == "Dryden"
        # Correlated PSDs (Van Steveren [PhD thesis] - Analysis of Aircraft Responses to Atmospheric Turbulence)
        Φuu = [π*σ^2/(1+aLΩx[i]^2+aLΩy[j]^2)^(5/2) * (1+aLΩx[i]^2+4*aLΩy[j]^2) for i in eachindex(aLΩx), j in eachindex(aLΩy)]
        Φuv = [π*σ^2/(1+aLΩx[i]^2+aLΩy[j]^2)^(5/2) * (-3*aLΩx[i]*aLΩy[j]/L) for i in eachindex(aLΩx), j in eachindex(aLΩy)]
        Φvv = [π*σ^2/(1+aLΩx[i]^2+aLΩy[j]^2)^(5/2) * (1+4*aLΩx[i]^2+aLΩy[j]^2) for i in eachindex(aLΩx), j in eachindex(aLΩy)]
        Φww = [π*σ^2/(1+aLΩx[i]^2+aLΩy[j]^2)^(5/2) * (3*(aLΩx[i]^2+aLΩy[j]^2)) for i in eachindex(aLΩx), j in eachindex(aLΩy)]
        # Transfer functions                                               
        H11 = @. sqrt(Φuu)
        H21 = @. Φuv/H11
        H22 = @. sqrt(Φvv-Φuv^2/Φuu)
        H33 = @. sqrt(Φww)
    end

    # Set random seed
    Random.seed!(seed)

    # White noise signals
    un = randn(Nx,Ny)
    vn = randn(Nx,Ny)
    wn = randn(Nx,Ny)

    # Gust velocity components in the spatial frequency domain
    if spectrum == "vK"
        Un = sqrt.(Φuu) .* fft(un)
        Vn = sqrt.(Φvv) .* fft(vn)
        Wn = sqrt.(Φww) .* fft(wn)
    elseif spectrum == "Dryden"        
        Un = H11 * fft(un)
        Vn = H21 * fft(un) + H22 * fft(vn)
        Wn = H33 * fft(wn)
    end

    # Gust velocity components in the time domain
    ug = real(ifft(Un))
    vg = real(ifft(Vn))
    wg = real(ifft(Wn))

    # Re-scale
    ug = ug*σ/rms(ug)
    vg = vg*σ/rms(vg)
    wg = wg*σ/rms(wg)

    # Spatial coordinates
    X = c0[1] .+ LinRange(0,gustWidth,Nx)
    Y = c0[2] .+ LinRange(0,gustLength,Ny)

    # Mapping functions
    itpU = Interpolate((X,Y), ug)
    itpV = Interpolate((X,Y), vg)
    itpW = Interpolate((X,Y), wg)

    # Gust velocity components as functions of spatial coordinate
    U = x::AbstractVector{<:Real} -> itpU(I2g12(x))
    V = x::AbstractVector{<:Real} -> itpV(I2g12(x))
    W = x::AbstractVector{<:Real} -> itpW(I2g12(x))

    # Set vector of gust velocity (in the inertial frame) as a function of position
    UGustInertial(x::AbstractVector{<:Real},t::Real) = isInside(I2g(x)) ? RT*[U(x); V(x); W(x)] : zeros(3)

    return Continuous2DSpaceGust(spectrum=spectrum,gustLength=gustLength,gustWidth=gustWidth,Nx=Nx,Ny=Ny,L=L,σ=σ,c0=c0,p=p,seed=seed,UGustInertial=UGustInertial,U=U,V=V,W=W)

end
export create_Continuous2DSpaceGust


# Computes a stochastic gust velocity array by coloring a white noise signal with the appropriate spectrum
function stochastic_gust_velocity_from_white_noise(spectrum::String,t::Vector{<:Real},URef::Real,L::Real=762,σ::Real=1)

    @assert spectrum in ["vK"]
    @assert maximum(abs.(diff(diff(t)))) < 1e-6 "t must be an evenly spaced vector"
    @assert L > 0
    @assert σ >= 0

    # Frequency vector
    # --------------------------------------------------------------------------
    # Size of time array
    N = length(t)
    # Time step
    Δt = t[2] - t[1]
    # Maximum detectable frequency (Nyquist frequency) [Hz]
    fmax = 1/(2*Δt)
    # Frequency vector [Hz]
    f = LinRange(0,fmax,N)

    # Transfer function and power spectra
    # --------------------------------------------------------------------------
    if spectrum == "vK"
        a = 1.339
        T = a*L/URef
        C = σ * sqrt(T/(a*π))
        Gu = @. C * (sqrt(2)/(1+1*im*f*T)^(5/6))
        Gvw = @. C * ((1+2*sqrt(2/3)*1*im*f*T)/(1+1*im*f*T)^(11/6))
    end
    Φu = @. abs(Gu)^2
    Φv = Φw = @. abs(Gvw)^2

    # FFT and PSD
    # --------------------------------------------------------------------------
    # White noise signals
    wn_u = randn(N)
    wn_v = randn(N)
    wn_w = randn(N)
    # FFT of white noise signals
    wnFFT_u = fft(wn_u)
    wnFFT_v = fft(wn_v)
    wnFFT_w = fft(wn_w)
    # FFT of the gust velocity components ("color" white noise with the power spectra)
    UFFT = @. wnFFT_u * sqrt(Φu)
    VFFT = @. wnFFT_v * sqrt(Φv)
    WFFT = @. wnFFT_w * sqrt(Φw)

    # Gust velocity components in time domain (re-scale to have the exact input RMS)
    # --------------------------------------------------------------------------
    U = real(ifft(UFFT))
    V = real(ifft(VFFT))
    W = real(ifft(WFFT))
    U = U*σ/rms(U)
    V = V*σ/rms(V)
    W = W*σ/rms(W)

    return U,V,W,f,Φu,Φv,Φw
end