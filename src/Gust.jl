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
- `initialTime::Real` = time when the gust begins
- `duration::Real` = duration of the gust
- `convectiveVelocity::Real` = convective velocity of the gust (in longitudinal direction)
- `verticalVelocity::Real` = peak "vertical" velocity of the gust (in lift direction)
- `p::Vector{<:Real}` = Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
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
    UGustInertial = t -> initialTime<t<finalTime ? RT*[0; convectiveVelocity; verticalVelocity] : zeros(3)

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
- `initialTime::Real` = time when the gust begins
- `duration::Real` = duration of the gust
- `convectiveVelocity::Real` = convective velocity of the gust (in longitudinal direction)
- `verticalVelocity::Real` = peak "vertical" velocity of the gust (in lift direction)
- `p::Vector{<:Real}` = Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
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
    UGustInertial = t -> initialTime<t<finalTime ? RT*[0; convectiveVelocity; 1/2*verticalVelocity*(1-cos(2π*(t-initialTime)/duration))] : zeros(3)

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
    L::Real
    ωmin::Real
    ωmax::Real
    Uref::Real
    convectiveVelocity::Real
    σ::Real
    p::Vector{<:Real}
    seed::Int64
    plotPSD::Bool
    
    # Secondary (fixed or outputs of gust creation)
    isDefinedOverTime::Bool = true
    UGustInertial::Function
    finalTime::Real
    V::Function

end


"""
    create_Continuous1DGust(; kwargs...)

Creates a continuous 1D gust

# Keyword arguments
- `spectrum::String` = spectrum of the PSD: von Kármán ("vK") or Dryden ("Dryden")
- `generationMethod::String` = method for generation of gust velocity ("sinusoids" or "whiteNoise")
- `initialTime::Real` = time when the gust begins
- `duration::Real` = duration of the gust
- `generationDuration::Real` = duration of the gust considered for its generation
- `L::Real` = turbulence length scale [m]
- `ωmin::Real` = minimum frequency of the PSD [rad/s]
- `ωmax::Real` = maximum frequency of the PSD [rad/s]
- `Uref::Real` = reference airspeed
- `convectiveVelocity::Real` = convective velocity of the gust (in longitudinal direction)
- `σ::Real` = turbulence intensity (RMS) of "vertical" gust velocity component
- `p::Vector{<:Real}` = Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
- `seed::Int64` = seed for random numbers generation (reproducibility)
- `plotPSD::Bool` = flag to plot the PSD
"""
function create_Continuous1DGust(; spectrum::String="vK",generationMethod::String="sinusoids",initialTime::Real=0,duration::Real,generationDuration::Real=duration,L::Real=762,ωmin::Real=0,ωmax::Real=200π,Uref::Real,convectiveVelocity::Real=0,σ::Real,p::Vector{<:Real}=zeros(3),seed::Int64=123456,plotPSD::Bool=false)

    # Validate
    @assert spectrum in ["vK","Dryden"]
    @assert generationMethod in ["sinusoids","whiteNoise"]
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

    # Get "vertical" velocity according to generation method
    if generationMethod == "sinusoids"
        # Frequency resolution [rad/s]: this ensures a stable value on the RMS of the gust velocity profile
        dω = min(2π/generationDuration, (2π)^2/ωmax) 
        # Frequency vector [rad/s]
        ω = collect(ωmin:dω:ωmax)
        # Spatial frequency vector and resolution [rad/m]
        Ω = ω/Uref
        dΩ = dω/Uref    
        # PSD [(m/s)^2/(rad/s)]
        if spectrum == "vK"
            Φ = @. σ^2*L/π * ( (1+8/3*(1.339*L*Ω)^2) / (1+(1.339*L*Ω)^2)^(11/6) )
        elseif spectrum == "Dryden"
            Φ = @. σ^2*L/π * ( (1+3*(L*Ω)^2) / (1+(L*Ω)^2)^2 )
        end
        sqrt2dΩΦ = sqrt.(2*dΩ*Φ)
        # Random phase vector
        ψ = 2π*rand(length(ω))
        # "Vertical" gust velocity over time
        V = t -> sum( sqrt2dΩΦ .* cos.(ω*t.+ψ) )
        # Plot PSD, if applicable
        if plotPSD
            Δt = π/ωmax
            time = collect(0:Δt:generationDuration-Δt)
            f,_,yPSD = get_FFT_and_PSD(time,V.(time))
            plt = plot(xscale=:log10,yscale=:log10,xlims=[L*Ω[2],L*Ω[end]],xlabel="Normalized frequency (\$\\Omega L\$)",ylabel="Normalized PSD")
            scatter!(L*Ω,Φ/(σ^2*L/π),c=:green,ms=5,msw=0,label="definition")
            plot!(L*Ω,yPSD/(2π)/(σ^2*L/(π*Uref)),c=:black,lw=2,label="generated")
            display(plt)
        end
    elseif generationMethod == "whiteNoise"
        # Time vector 
        Δt = π/ωmax
        time = collect(0:Δt:generationDuration-Δt)
        # Gust velocity array
        Vwn,f,Φ = stochastic_gust_velocity_from_white_noise(spectrum,time,Uref,L,σ)
        # Mapping function
        itp = Interpolate(time, Vwn)
        # "Vertical" gust velocity over time
        V = t -> itp(t)
        # Plot PSD, if applicable
        if plotPSD
            f2,_,yPSD = get_FFT_and_PSD(time,Vwn)
            plt = plot(xscale=:log10,yscale=:log10,xlims=[f[2],f[end]],xlabel="frequency [Hz]",ylabel="Normalized PSD")
            plot!(f2,yPSD/(σ^2*L/(π*Uref)),lw=2,label="generated")
            plot!(f,Φ/(σ^2*L/(π*Uref)),lw=2,label="definition")
            display(plt)
        end
    end

    # Rotation tensor from the inertial basis to the gust basis
    R = rotation_tensor_E321(p)
    RT = R'

    # Final time of the gust
    finalTime = initialTime + duration

    # Set vector of gust velocity (in the inertial frame) function over time
    UGustInertial = t -> initialTime<t<finalTime ? RT*[0; convectiveVelocity; V(t)] : zeros(3)

    return Continuous1DGust(spectrum=spectrum,generationMethod=generationMethod,initialTime=initialTime,duration=duration,generationDuration=generationDuration,L=L,ωmin=ωmin,ωmax=ωmax,Uref=Uref,convectiveVelocity=convectiveVelocity,σ=σ,p=p,seed=seed,plotPSD=plotPSD,finalTime=finalTime,UGustInertial=UGustInertial,V=V)

end
export create_Continuous1DGust


#
# @with_kw mutable struct DiscreteSpaceGust <: Gust

#     DiscreteSpaceGust composite type

#
@with_kw mutable struct DiscreteSpaceGust <: Gust

    # Primary (necessary for gust creation)
    type::String
    length::Real
    width::Real
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
- `type::String` = type of gust
- `length::Real` = length of the gust (in longitudinal direction)
- `width::Real` = width of the gust (in spanwise direction)
- `convectiveVelocity::Real` = convective velocity of the gust (in longitudinal direction)
- `verticalVelocity::Real` = peak "vertical" velocity of the gust (in lift direction)
- `c0::Vector{<:Real}` = position vector of the front of the gust, resolved in the inertial basis
- `p::Vector{<:Real}` = Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
"""
function create_DiscreteSpaceGust(; type::String,length::Real,width::Real,convectiveVelocity::Real=0,verticalVelocity::Real,c0::Vector{<:Real}=zeros(3),p::Vector{<:Real}=zeros(3))

    # Validate
    @assert type in ["SharpEdged","OneMinusCosine","DARPA"]
    @assert length > 0
    @assert width > 0
    @assert size(c0) == (3,)
    @assert size(p) == (3,)

    # Rotation tensor from the inertial basis to the gust basis
    R = rotation_tensor_E321(p)
    RT = R'

    # Normalized width coordinate (ranges from -1/2 to 1/2 for -width/2 <= x[1] <= width/2)
    w = x -> (x[1]-c0[1])/width

    # Normalized length coordinate (ranges from 0 to 1 for 0 <= x[2] <= length)
    l = x -> (x[2]-c0[2])/length

    # Gust front longitudinal coordinate as a function of time
    c(t) = c0[2] + convectiveVelocity*t

    # TF function for a position vector being inside the gust at a given time
    isInside(x,t) = -width/2 <= (x[1]-c0[1]) <= width/2 && 0 <= x[2]-c(t) <= length

    # Vector transformation from inertial to gust basis
    I2g = x -> RT * x
    
    # Set "vertical" velocity function according to gust type
    if type == "SharpEdged"
        V = x -> verticalVelocity
    elseif type == "OneMinusCosine"
        V = x -> 1/2*verticalVelocity*(1-cos(2π*l(I2g(x))))
    elseif type == "DARPA"
        V = x -> -1/2*verticalVelocity*cos(2π*w(I2g(x)))*(1-cos(2π*l(I2g(x))))
    end

    # Set vector of gust velocity (in the inertial frame) as a function of position and time
    UGustInertial(x,t) = isInside(I2g(x),t) ? RT*[0; convectiveVelocity; V(x)] : zeros(3)

    return DiscreteSpaceGust(type=type,length=length,width=width,convectiveVelocity=convectiveVelocity,verticalVelocity=verticalVelocity,c0=c0,p=p,UGustInertial=UGustInertial)

end
export create_DiscreteSpaceGust


#
# @with_kw mutable struct Continuous1DSpaceGust <: Gust

#     Continuous1DSpaceGust composite type

#
@with_kw mutable struct Continuous1DSpaceGust <: Gust

    # Primary (necessary for gust creation)
    spectrum::String
    length::Real
    N::Int64
    L::Real
    σ::Real
    c0::Vector{<:Real}
    p::Vector{<:Real}
    seed::Int64
    
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
- `spectrum::String` = spectrum of the PSD: von Kármán ("vK") or Dryden ("Dryden")
- `length::Real` = length of the gust (in longitudinal direction)
- `N::Int64` = number of nodes for length discretization
- `L::Real` = turbulence length scale
- `σ::Real` = turbulence intensity (RMS) of "vertical" gust velocity component
- `c0::Vector{<:Real}` = position vector of the front of the gust, resolved in the inertial basis
- `p::Vector{<:Real}` = Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
- `seed::Int64` = seed for random numbers generation (reproducibility)
"""
function create_Continuous1DSpaceGust(; spectrum::String,length::Real,N::Int64=1001,L::Real=762,σ::Real=1,c0::Vector{<:Real}=zeros(3),p::Vector{<:Real}=zeros(3),seed::Int64=123456)

    # Validate
    @assert spectrum in ["vK","Dryden"]
    @assert length > 0
    @assert N > 0
    @assert L > 0
    @assert size(c0) == (3,)
    @assert size(p) == (3,)

    # Rotation tensor from the inertial basis to the gust basis
    R = rotation_tensor_E321(p)
    RT = R'

    # TF function for a position vector being inside the gust
    isInside(x) = 0 <= x[2]-c0[2] <= length

    # Vector transformation from inertial to gust basis
    I2g = x -> RT * x

    # Second (longitudinal) component of the above
    I2g2 = x -> I2g(x)[2]

    # Non-dimensional spatial frequency vector
    LΩ = 2π*L/length * LinRange(0,div(N-1,2),N)

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

    # Gust velocity components in the time domain
    ug = real(ifft(Un))
    vg = real(ifft(Vn))
    wg = real(ifft(Wn))

    # Re-scale
    ug = ug*σ/rms(ug)
    vg = vg*σ/rms(vg)
    wg = wg*σ/rms(wg)

    # Spatial coordinate 
    X = c0[2] .+ LinRange(0,length,N)

    # Mapping functions
    itpU = Interpolate(X, ug)
    itpV = Interpolate(X, vg)
    itpW = Interpolate(X, wg)

    # Gust velocity components as functions of spatial coordinate
    U = x -> itpU(I2g2(x))
    V = x -> itpV(I2g2(x))
    W = x -> itpW(I2g2(x))

    # Set vector of gust velocity (in the inertial frame) as a function of position
    UGustInertial(x,t) = isInside(I2g(x)) ? RT*[U(x); V(x); W(x)] : zeros(3)

    return Continuous1DSpaceGust(spectrum=spectrum,length=length,N=N,L=L,σ=σ,c0=c0,p=p,seed=seed,UGustInertial=UGustInertial,U=U,V=V,W=W)

end
export create_Continuous1DSpaceGust


#
# @with_kw mutable struct Continuous2DSpaceGust <: Gust

#     Continuous2DSpaceGust composite type

#
@with_kw mutable struct Continuous2DSpaceGust <: Gust

    # Primary (necessary for gust creation)
    spectrum::String
    length::Real
    width::Real
    Nx::Int64
    Ny::Int64
    L::Real
    σ::Real
    c0::Vector{<:Real}
    p::Vector{<:Real}
    seed::Int64
    
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
- `spectrum::String` = spectrum of the PSD: von Kármán ("vK") or Dryden ("Dryden")
- `length::Real` = length of the gust (in longitudinal direction)
- `width::Real` = width of the gust (in lateral direction)
- `Nx::Int64` = number of nodes for length discretization
- `Ny::Int64` = number of nodes for width discretization
- `L::Real` = turbulence length scale
- `σ::Real` = turbulence intensity (RMS) of "vertical" gust velocity component
- `c0::Vector{<:Real}` = position vector of the front of the gust, resolved in the inertial basis
- `p::Vector{<:Real}` = Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
- `seed::Int64` = seed for random numbers generation (reproducibility)
"""
function create_Continuous2DSpaceGust(; spectrum::String,length::Real,width::Real,Nx::Int64=101,Ny::Int64=101,L::Real=762,σ::Real=1,c0::Vector{<:Real}=zeros(3),p::Vector{<:Real}=zeros(3),seed::Int64=123456)

    # Validate
    @assert spectrum in ["vK","Dryden"]
    @assert length > 0
    @assert width > 0
    @assert Nx > 0
    @assert Ny > 0
    @assert L > 0
    @assert size(c0) == (3,)
    @assert size(p) == (3,)

    # Rotation tensor from the inertial basis to the gust basis
    R = rotation_tensor_E321(p)
    RT = R'

    # TF function for a position vector being inside the gust
    isInside(x) = -width/2 <= x[1]-c0[1] <= width/2 && 0 <= x[2]-c0[2] <= length

    # Vector transformation from inertial to gust basis
    I2g = x -> RT * x

    # Longitudinal and lateral components of the above
    I2g12 = x -> I2g(x)[1:2]

    # Constant
    if spectrum == "vK"
        a = 1.339
    elseif spectrum == "Dryden"
        a = 1
    end

    # Non-dimensional spatial frequency vectors
    aLΩx = 2π*a*L/width * LinRange(0,div(Nx-1,2),Nx)
    aLΩy = 2π*a*L/length * LinRange(0,div(Ny-1,2),Ny)

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
    X = c0[1] .+ LinRange(0,width,Nx)
    Y = c0[2] .+ LinRange(0,length,Ny)

    # Mapping functions
    itpU = Interpolate((X,Y), ug)
    itpV = Interpolate((X,Y), vg)
    itpW = Interpolate((X,Y), wg)

    # Gust velocity components as functions of spatial coordinate
    U = x -> itpU(I2g12(x))
    V = x -> itpV(I2g12(x))
    W = x -> itpW(I2g12(x))

    # Set vector of gust velocity (in the inertial frame) as a function of position
    UGustInertial(x,t) = isInside(I2g(x)) ? RT*[U(x); V(x); W(x)] : zeros(3)

    return Continuous2DSpaceGust(spectrum=spectrum,length=length,width=width,Nx=Nx,Ny=Ny,L=L,σ=σ,c0=c0,p=p,seed=seed,UGustInertial=UGustInertial,U=U,V=V,W=W)

end
export create_Continuous2DSpaceGust


# Computes a stochastic gust velocity array by coloring a white noise signal with the appropriate spectrum
function stochastic_gust_velocity_from_white_noise(spectrum::String,t::Vector{<:Real},U::Real,L::Real=762,σ::Real=1)

    @assert spectrum in ["vK"]
    @assert maximum(abs.(diff(diff(t)))) < 1e-10 "t must be an evenly spaced vector"
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

    # Transfer function and power spectrum
    # --------------------------------------------------------------------------
    if spectrum == "vK"
        a = 1.339
        T = a*L/U
        C = σ * sqrt(T/(a*π))
        G = @. C * ((1+2*sqrt(2/3)*1*im*f*T)/(1+1*im*f*T)^(11/6))
    end
    Φ = @. abs(G)^2

    # FFT and PSD
    # --------------------------------------------------------------------------
    # White noise signal
    wn = randn(N)
    # FFT of white noise
    wnFFT = fft(wn)
    # FFT of the gust velocity ("color" white noise with the power spectrum)
    VFFT = @. wnFFT * sqrt(Φ)

    # Gust velocity in time domain (re-scale to have the exact input RMS)
    # --------------------------------------------------------------------------
    V = real(ifft(VFFT))
    V = V*σ/rms(V)

    return V,f,Φ
end