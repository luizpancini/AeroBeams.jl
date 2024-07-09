using SparseArrays, LinearAlgebra, BenchmarkTools, Arpack

N = 1000
A::SparseMatrixCSC{Float64,Int64} = sprandn(N,N,0.5)
B::SparseMatrixCSC{Float64,Int64} = sprandn(N,N,0.5)
B2 = Matrix(B)

function fun1(A,B)
    return eigen(A\Matrix(B))
end

function fun2(A,B)
    return eigen(A\B)
end

function fun3(A,B)
    return eigen(A\Array(B))
end

@btime fun1(A,B)
@btime fun2(A,B2)
@btime fun3(A,B)