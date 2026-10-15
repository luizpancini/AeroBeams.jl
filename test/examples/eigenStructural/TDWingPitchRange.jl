using AeroBeams

# Pitch angle range [rad]
θRange = LinRange(0,π/2,21)

# Discretization
nElem = 20

# Initialize model
TDWingPitchRange = create_TDWing(nElem=nElem)

# Number of modes to track
nModes = 4

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :arnoldi
nev = 40

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(θRange))
tip_u2 = Array{Float64}(undef,length(θRange))
tip_u3 = Array{Float64}(undef,length(θRange))
tip_twist = Array{Float64}(undef,length(θRange))

# Loop configurations
for (i,θ) in enumerate(θRange)
    display("Solving for root angle = $(round(θ*180/π,digits=3)) deg")
    # Update wing pitch angle
    update_beam_orientation!(TDWingPitchRange; beam=TDWingPitchRange.beams[1],p0=transform_to_WM([0; 0; θ],from=:E321))
    # Create and solve problem
    problem = create_EigenProblem(model=TDWingPitchRange,eigenmethod=eigenmethod,nev=nev)
    solve!(problem)
    # Eigensolution
    eigensols[i] = problem.eigensolution
    # Get outputs
    lastNodeFinalState = problem.loadstepHistory.states[end][nElem].nodal
    tip_u2[i] = lastNodeFinalState.u2_b[2]
    tip_u3[i] = lastNodeFinalState.u2_b[3]
    tip_twist[i] = twist_angle(lastNodeFinalState.p2_b)
end

# Apply mode tracking
freqs,_ = mode_tracking(eigensols)

# Get frequencies of oscillatory modes in Hz
freqs,_ = filter_tracked_modes(rads2Hz.(freqs); fmin=eps(), maxTracks=nModes)

# Load reference solutions
u2_exp = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "u2_exp.txt"))
u3_exp = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "u3_exp.txt"))
th_exp = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "th_exp.txt"))
u2_num = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "u2_num.txt"))
u3_num = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "u3_num.txt"))
th_num = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "th_num.txt"))
freqs_exp = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "freqs_exp.txt"))
freq1_num = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "freq1_num.txt"))
freq2_num = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "freq2_num.txt"))
freq4_num = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "freq4_num.txt"))
freqs_num = [freq1_num, freq2_num, freq4_num]

## Plots
# Set paths
relPath = "/test/outputs/figures/TDWingPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:darkrainbow, 3, categorical=true)
lw = 2
ms = 5
msw = 0
gr()

# Tip flapwise displacement
plt1 = plot(xlabel="Root angle [deg]", ylabel="Tip flapwise displacement [m]", xticks=collect(0:15:90))
plot!(rad2deg.(θRange), -tip_u3, c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!(u3_num[1,:], u3_num[2,:], c=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!(u3_exp[1,:], u3_exp[2,:], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
display(plt1)
savefig(string(absPath,"/TDWingPitchRange_u3.pdf"))

# Tip chordwise displacement
plt2 = plot(xlabel="Root angle [deg]", ylabel="Tip chordwise displacement [m]", xticks=collect(0:15:90))
plot!(rad2deg.(θRange),-tip_u2, c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!(u2_num[1,:],u2_num[2,:], c=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!(u2_exp[1,:],u2_exp[2,:], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
display(plt2)
savefig(string(absPath,"/TDWingPitchRange_u2.pdf"))

# Tip twist
plt3 = plot(xlabel="Root angle [deg]", ylabel="Tip twist [deg]", xticks=collect(0:15:90), legend=:bottom)
plot!(rad2deg.(θRange),-tip_twist, c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!(th_num[1,:],th_num[2,:], c=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!(th_exp[1,:],th_exp[2,:], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
display(plt3)
savefig(string(absPath,"/TDWingPitchRange_twist.pdf"))

# Structural frequencies 
plt4 = plot(xlabel="Root angle [deg]", ylabel="Frequency [Hz]", legend=:outertop, xticks=collect(0:15:90), ylims=[0,30])
modeLabels = ["Flapwise bending" "Chordwise bending" "Torsion"]
plot!([NaN],[NaN], c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN],[NaN], c=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
for (m,mode) in enumerate([1,2,4])
    plot!(rad2deg.(θRange), freqs[mode,:], c=colors[m], lw=lw, ls=:solid, label=false)
    plot!(freqs_num[m][1,:],freqs_num[m][2,:], c=colors[m], lw=lw, ls=:dash, label=false)
    scatter!(freqs_exp[1,:],freqs_exp[1+m,:], mc=colors[m], ms=ms, msw=msw, label=false)
    isfinite(freqs[mode,1]) && annotate!(30, freqs[mode,1], text(modeLabels[m], :bottom, colors[m]))
end
display(plt4)
savefig(string(absPath,"/TDWingPitchRange_freqs.pdf"))

println("Finished TDWingPitchRange.jl")
