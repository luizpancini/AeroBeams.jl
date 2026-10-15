using AeroBeams

# Sweep angle range [rad]
ΛRange = deg2rad.(0:5:30)

# Flag for ad hoc sectional stiffness corrections with sweep angle
sweepStructuralCorrections = true

# Gravity vector
g = zeros(3)

# Flag for upright configuration
upright = false

# Number of modes to track
nModes = 9

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :fulldecomp
nev = 50

# Initialize outputs
problem = Array{EigenProblem}(undef,length(ΛRange))
eigensols = Array{EigenSolution}(undef,length(ΛRange))

# Loop sweep angle
for (i,Λ) in enumerate(ΛRange)
    display("Solving for Λ = $(round(Int,rad2deg(Λ))) deg")
    # Model
    sweptPazyModal,_ = create_Pazy(sweepAngle=Λ,upright=upright,gravityVector=g,sweepStructuralCorrections=sweepStructuralCorrections)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=sweptPazyModal,eigenmethod=eigenmethod,nev=nev)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
end

# Apply mode tracking
freqs,_ = mode_tracking(eigensols)

# Get frequencies of oscillatory modes in Hz
freqs,_ = filter_tracked_modes(rads2Hz.(freqs); fmin=eps(), maxTracks=nModes)

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptPazyModal"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot 10 deg swept mode shapes
i2plot = 3
plt_swept_modeShapes = plot_mode_shapes(problem[i2plot],nModes=nModes,scale=0.2,view=(30,30),legendPos=(-0.25,0.99),save=true,savePath=string(relPath,"/PazySweepAngleRangeModal_modeShapes.pdf"))
display(plt_swept_modeShapes)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 4
msw = 0
gr()

# Frequencies vs sweep angle
plt1 = plot(xlabel="Sweep angle [deg]", ylabel="Frequency [Hz]", ylims=[0,120])
for mode in 1:nModes
    plot!(rad2deg.(ΛRange), freqs[mode,:], c=colors[mode], lw=lw, marker=:circle, ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/sweptPazyModal.pdf"))

println("Finished sweptPazyModal.jl")