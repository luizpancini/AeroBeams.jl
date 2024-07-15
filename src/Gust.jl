abstract type Gust end

"""
@with_kw mutable struct SharpEdged <: Gust

    SharpEdged composite type

# Fields
- `initialTime::Number` = time when the gust begins
- `duration::Number` = duration of the gust
- `convectiveVelocity::Number` = convective velocity of the gust (in longitudinal direction)
- `verticalVelocity::Number` = peak "vertical" velocity of the gust (in lift direction)
- `p::Vector{<:Number}` = Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
- `isDefinedOverTime::Bool` = TF for being a gust defined over time (true)
- `UGustInertial::Function` = inertial gust velocity vector, as a function of time
- `finalTime::Number` = time when the gust ends
"""
@with_kw mutable struct SharpEdged <: Gust

    # Primary (necessary for gust creation)
    initialTime::Number
    duration::Number
    convectiveVelocity::Number
    verticalVelocity::Number
    p::Vector{<:Number}
    
    # Secondary (fixed or outputs of gust creation)
    isDefinedOverTime::Bool = true
    UGustInertial::Function
    finalTime::Number

end

# Constructor
function create_SharpEdgedGust(; initialTime::Number=0,duration::Number=Inf,convectiveVelocity::Number=0,verticalVelocity::Number,p::Vector{<:Number}=zeros(3))

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
    UGustInertial = t -> ifelse(initialTime<t<finalTime, RT*[0; convectiveVelocity; verticalVelocity], zeros(3))

    return SharpEdged(initialTime=initialTime,duration=duration,convectiveVelocity=convectiveVelocity,verticalVelocity=verticalVelocity,p=p,finalTime=finalTime,UGustInertial=UGustInertial)

end
export create_SharpEdgedGust


"""
@with_kw mutable struct OneMinusCosine <: Gust

    OneMinusCosine composite type

# Fields
- `initialTime::Number` = time when the gust begins
- `duration::Number` = duration of the gust
- `convectiveVelocity::Number` = convective velocity of the gust (in longitudinal direction)
- `verticalVelocity::Number` = peak "vertical" velocity of the gust (in lift direction)
- `p::Vector{<:Number}` = Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
- `isDefinedOverTime::Bool` = TF for being a gust defined over time (true)
- `UGustInertial::Function` = inertial gust velocity vector, as a function of time
- `finalTime::Number` = time when the gust ends
"""
@with_kw mutable struct OneMinusCosine <: Gust

    # Primary (necessary for gust creation)
    initialTime::Number
    duration::Number
    convectiveVelocity::Number
    verticalVelocity::Number
    p::Vector{<:Number}
    
    # Secondary (fixed or outputs of gust creation)
    isDefinedOverTime::Bool = true
    UGustInertial::Function
    finalTime::Number

end

# Constructor
function create_OneMinusCosineGust(; initialTime::Number=0,duration::Number,convectiveVelocity::Number=0,verticalVelocity::Number,p::Vector{<:Number}=zeros(3))

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
    UGustInertial = t -> ifelse(initialTime<t<finalTime, RT*[0; convectiveVelocity; 1/2*verticalVelocity*(1-cos(2π*(t-initialTime)/duration))], zeros(3))

    return OneMinusCosine(initialTime=initialTime,duration=duration,convectiveVelocity=convectiveVelocity,verticalVelocity=verticalVelocity,p=p,finalTime=finalTime,UGustInertial=UGustInertial)

end
export create_OneMinusCosineGust


"""
@with_kw mutable struct DARPA <: Gust

    DARPA composite type

# Fields
- `length::Number` = length of the gust (in longitudinal direction)
- `width::Number` = width of the gust (in spanwise direction)
- `convectiveVelocity::Number` = convective velocity of the gust (in longitudinal direction)
- `verticalVelocity::Number` = peak "vertical" velocity of the gust (in lift direction)
- `c0::Vector{<:Number}` = position vector of the front of the gust, resolved in the inertial basis
- `p::Vector{<:Number}` = Euler rotation parameters (3-2-1 sequence) from the inertial basis to the gust basis
- `isDefinedOverTime::Bool` = TF for being a gust defined over time (false)
- `UGustInertial::Function` = inertial gust velocity vector, as a function of position
"""
@with_kw mutable struct DARPA <: Gust

    # Primary (necessary for gust creation)
    length::Number
    width::Number
    convectiveVelocity::Number
    verticalVelocity::Number
    c0::Vector{<:Number}
    p::Vector{<:Number}
    
    # Secondary (fixed or outputs of gust creation)
    isDefinedOverTime::Bool = false
    UGustInertial::Function

end

# Constructor
function create_DARPAGust(; length::Number,width::Number,convectiveVelocity::Number=0,verticalVelocity::Number,c0::Vector{<:Number}=zeros(3),p::Vector{<:Number}=zeros(3))

    # Validate
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
    
    # Set vector of gust velocity (in the inertial frame) as a function of position
    UGustInertial(x,t) = ifelse(isInside(I2g(x),t), RT*[0; convectiveVelocity; -1/2*verticalVelocity*cos(2π*w(I2g(x)))*(1-cos(2π*l(I2g(x))))], zeros(3))

    return DARPA(length=length,width=width,convectiveVelocity=convectiveVelocity,verticalVelocity=verticalVelocity,c0=c0,p=p,UGustInertial=UGustInertial)

end
export create_DARPAGust