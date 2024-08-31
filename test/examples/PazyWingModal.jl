using AeroBeams

# Number of modes
nModes = 5

# Initialize outputs
problem = Array{EigenProblem}(undef,2)
freqs = Array{Vector{Float64}}(undef,2)

# Loop position configurations
for (i,upright) in enumerate([false,true])
    # Model
    PazyWingModal,_ = create_Pazy(upright=upright)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=PazyWingModal,nModes=nModes)
    solve!(problem[i])
    # Get frequencies in Hz
    freqs[i] = problem[i].frequenciesOscillatory/(2π)
end

# Load reference data
modalFreqsGVT = [4.39 29.80 41.00 82.50 NaN; 4.26 28.50 42.00 81.50 60.70]
modalFreqsUMNAST = [4.19 28.49 41.58 83.03 100.07; 4.19 28.49 41.88 83.06 105.89]

# Compare solutions
ϵ_rel_horz_GVT = freqs[1] ./ modalFreqsGVT[1,:] .- 1
ϵ_rel_vert_GVT = freqs[2] ./ modalFreqsGVT[2,:] .- 1
ϵ_rel_horz_UMNAST = freqs[1] ./ modalFreqsUMNAST[1,:] .- 1
ϵ_rel_vert_UMNAST = freqs[2] ./ modalFreqsUMNAST[2,:] .- 1

println("Relative erros for horizontal wing:")
println("GVT: $ϵ_rel_horz_GVT")
println("UMNAST: $ϵ_rel_horz_UMNAST")
println("Relative erros for vertical wing:")
println("GVT: $ϵ_rel_vert_GVT")
println("UMNAST: $ϵ_rel_vert_UMNAST")

println("Finished PazyWingModal.jl")