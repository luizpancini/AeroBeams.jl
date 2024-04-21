using Parameters, LinearAlgebra

@with_kw mutable struct MyStruct
    C::Vector{Matrix{<:Number}}
end

function create_mystruct(;C::Vector{<:Matrix{<:Number}})
    return MyStruct(C)
end

c = diagm([1.0,1.0])

s = create_mystruct(C=[c])

s.C