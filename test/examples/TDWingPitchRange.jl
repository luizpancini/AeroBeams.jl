using AeroBeams, LinearAlgebra, DelimitedFiles

# Discretization
nElem = 20

# Loop configurations
θRange = LinRange(0,π/2,21) 
freqs = Array{Vector{Float64}}(undef,length(θRange))
tip_u2 = Array{Float64}(undef,length(θRange))
tip_u3 = Array{Float64}(undef,length(θRange))
tip_twist = Array{Float64}(undef,length(θRange))
for (i,θ) in enumerate(θRange)
    display("Solving for root angle = $(round(θ*180/π,digits=3)) deg")
    # Model
    TDWingPitchRange = create_TDWing(θ=θ,nElem=nElem)
    # Create and solve problem
    global problem = create_EigenProblem(model=TDWingPitchRange,nModes=4,frequencyFilterLimits=[0.1,Inf64],normalizeModeShapes=true)
    solve!(problem)
    # Get outputs
    freqs[i] = problem.frequenciesOscillatory/(2π)
    tip_u2[i] = problem.nodalStatesOverσ[end][nElem].u_n2_b[2]
    tip_u3[i] = problem.nodalStatesOverσ[end][nElem].u_n2_b[3]
    tip_p = problem.nodalStatesOverσ[end][nElem].p_n2_b
    R,_ = rotation_tensor_WM(tip_p)
    Δ = R*[0; 1; 0]
    tip_twist[i] = asind(Δ[3])
end

# Load reference solutions
u2_exp = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "u2_exp.txt"))
u3_exp = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "u3_exp.txt"))
th_exp = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "th_exp.txt"))
u2_num = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "u2_num.txt"))
u3_num = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "u3_num.txt"))
th_num = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "th_num.txt"))
freqs_exp = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "freqs_exp.txt"))
freq1_num = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "freq1_num.txt"))
freq2_num = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "freq2_num.txt"))
freq4_num = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "freq4_num.txt"))
freqs_num = Array{Matrix{Float64}}(undef,3)
freqs_num[1] = freq1_num
freqs_num[2] = freq2_num
freqs_num[3] = freq4_num

println("Finished TDWingPitchRange.jl")