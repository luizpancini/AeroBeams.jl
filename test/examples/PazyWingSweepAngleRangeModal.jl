# Not a currently working example: need to change the sectional stiffness and inertia matrices with the sweep angle

using AeroBeams

# Number of modes
nModes = 3

# Gravity
g = 0

# Sweep angle [rad] range
ΛRange = [0,10,30]*π/180

# Initialize outputs
problem = Array{EigenProblem}(undef,length(ΛRange))
freqs = Array{Vector{Float64}}(undef,length(ΛRange))

# Loop sweep angle
for (i,Λ) in enumerate(ΛRange)
    # Model
    PazyWingSweepAngleRangeModal,_ = create_Pazy(Λ=Λ,upright=false,g=g)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=PazyWingSweepAngleRangeModal,nModes=nModes)
    solve!(problem[i])
    # Get frequencies in Hz
    freqs[i] = problem[i].frequenciesOscillatory/(2π)
end

# Separate frequencies by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(ΛRange)]
end

using Plots, ColorSchemes

# Plot configurations
colors = get(colorschemes[:darkrainbow], LinRange(0, 1, nModes))
lw = 2
gr()

# Frequencies vs sweep angle
plt1 = plot(xlabel="Sweep angle [deg]", ylabel="Frequency [Hz]")
for mode in 1:nModes
    plot!(ΛRange, modeFrequencies[mode], c=colors[mode], lw=lw,  label=false)
end
display(plt1)
savefig(string(absPath,"/PazyWingSweepAngleRangeModal.pdf"))

println("Finished PazyWingSweepAngleRangeModal.jl")