using AeroBeams, Plots, ColorSchemes

# Number of modes
nModes = 5

# Gravity
g = 0

# Sweep angle [rad] range
ΛRange = vcat(0:5:30)*π/180

# Flag for ad hoc sectional stiffness corrections with sweep angle
sweepStructuralCorrections = true

# Initialize outputs
problem = Array{EigenProblem}(undef,length(ΛRange))
freqs = Array{Vector{Float64}}(undef,length(ΛRange))

# System solver
absTol = 1e-7
NR = create_NewtonRaphson(displayStatus=true,absoluteTolerance=absTol)

# Loop sweep angle
for (i,Λ) in enumerate(ΛRange)
    display("Solving for Λ = $(round(Λ*180/π)) deg")
    # Model
    sweptPazyModal,_ = create_Pazy(Λ=Λ,upright=false,g=g,sweepStructuralCorrections=sweepStructuralCorrections)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=sweptPazyModal,nModes=nModes,systemSolver=NR)
    solve!(problem[i])
    # Get frequencies in Hz
    freqs[i] = problem[i].frequenciesOscillatory/(2π)
end

# Separate frequencies by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(ΛRange)]
end

# Set paths
relPath = "/dev/sweptPazy/General/outputs/sweptPazyModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot swept configuration
plt_swept_modeShapes = plot_mode_shapes(problem[end],view=(30,30),save=true,savePath=string(relPath,"/PazySweepAngleRangeModal_modeShapes.pdf"))
display(plt_swept_modeShapes)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 4
msw = 0
gr()

# Frequencies vs sweep angle
plt1 = plot(xlabel="Sweep angle [deg]", ylabel="Frequency [Hz]", ylims=[0,50])
for mode in 1:nModes
    plot!(ΛRange*180/π, modeFrequencies[mode], c=colors[mode], lw=lw, marker=:circle, ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/sweptPazyModal.pdf"))

println("Finished sweptPazyModal.jl")