using BenchmarkTools, StaticArrays

function compute_all_at_once(f,time)
    return f.(time)
end

function compute_one_at_a_time(f,time)
    fun = SVector{length(time),Float64}
    for i in eachindex(time)
        fun[i] = f(time[i])
    end
    return fun
end

f = t -> exp.(t).*t.^2

time = collect(0:0.001:10)

@btime compute_all_at_once(f,time)
@btime compute_one_at_a_time(f,time)


