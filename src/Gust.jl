abstract type Gust end

"""
@with_kw mutable struct SharpEdged <: Gust

    SharpEdged composite type

# Fields
- 
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
    UGustInertial = t -> RT*[0; convectiveVelocity; ifelse(initialTime<t<finalTime, verticalVelocity, 0)]

    return OneMinusCosine(initialTime=initialTime,duration=duration,convectiveVelocity=convectiveVelocity,verticalVelocity=verticalVelocity,p=p,finalTime=finalTime,UGustInertial=UGustInertial)

end
export create_SharpEdgedGust


"""
@with_kw mutable struct OneMinusCosine <: Gust

    OneMinusCosine composite type

# Fields
- 
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
    UGustInertial = t -> RT*[0; convectiveVelocity; ifelse(initialTime<t<finalTime, 1/2*verticalVelocity*(1-cos(2π*(t-initialTime)/duration)), 0)]

    return OneMinusCosine(initialTime=initialTime,duration=duration,convectiveVelocity=convectiveVelocity,verticalVelocity=verticalVelocity,p=p,finalTime=finalTime,UGustInertial=UGustInertial)

end
export create_OneMinusCosineGust