using AeroBeams

# Sweep angle range [rad]
ΛRange = deg2rad.([0,10,20,30])

# Tip load [N]
F3TipRange = vcat(0:1:35)

# Flag for ad hoc sectional stiffness corrections with sweep angle
sweepStructuralCorrections = true

# Gravity
g = 0

# Number of modes to track
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :fulldecomp
nev = 40

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(ΛRange),length(F3TipRange))
freqsOsc = Array{Matrix{Float64}}(undef,length(ΛRange))
tipOOP = Array{Float64}(undef,length(ΛRange),length(F3TipRange))
tip_twist = Array{Float64}(undef,length(ΛRange),length(F3TipRange))

# Loop sweep angle
for (i,Λ) in enumerate(ΛRange)
    # Initialize guess solution
    x0 = zeros(0)
    for (j,F3) in enumerate(F3TipRange)
        display("Solving for Λ = $(round(Int,rad2deg(Λ))) deg, tip force = $F3 N")
        # Dummy beam for the BC
        dummyBeam = create_Beam(length=1,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1)])
        # Tip force BC
        tipForce = create_BC(name="tipForce",beam=dummyBeam,node=nElem+1,conditions=Dict(:F3b=>F3))
        # Model
        sweptPazyModalTipDispRange,_ = create_Pazy(sweepAngle=Λ,upright=false,gravityVector=[0,0,-g],sweepStructuralCorrections=sweepStructuralCorrections,additionalBCs=[tipForce])
        # Create and solve problem
        problem = create_EigenProblem(model=sweptPazyModalTipDispRange,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(problem)
        # Eigensolution
        eigensols[i,j] = problem.eigensolution
        # Displacement outputs
        tipNodeFinalState = problem.loadstepHistory.states[end][nElem].nodal
        tip_twist[i,j] = twist_angle(tipNodeFinalState.p2_b)
        tipOOP[i,j] = tipNodeFinalState.u2_b[3]
        # Set initial guess solution for next iteration
        x0 = problem.systemArrays.x
    end
    # Apply mode tracking
    freqs,_ = mode_tracking(eigensols[i,:])
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[i],_ = filter_tracked_modes(rads2Hz.(freqs); fmin=eps(),maxTracks=nModes)
end

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptPazyModalTipDispRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lfs = 12
lw = 2
ms = 5
ls = [:solid :dash :dot :dashdot]
gr()

# Frequencies vs sweep angle
plt1 = plot(xlabel="Normalized tip displacement", ylabel="Frequency [Hz]", xlims=[0,0.5], ylims=[0,50], legendfontsize=lfs, colorbar_title="Mode")
for (i,Λ) in enumerate(ΛRange)
    plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[i], label="\$\\Lambda=\$$(round.(Int,rad2deg(Λ)))\$^\\circ\$")
end
for (i,Λ) in enumerate(ΛRange)
    for mode in 1:nModes
        plot!(tipOOP[i,:]/L, freqsOsc[i][mode,:], c=colors, lw=lw, ls=ls[i], lz=mode, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/sweptPazyModalTipDispRange.pdf"))

# Tip displacement vs. sweep angle for largest load
plt_OOP = plot(xlabel="Sweep angle [deg]", ylabel="Tip \$u_3/L\$")
plot!(rad2deg.(ΛRange), tipOOP[:,end]/L, c=:black, lw=lw, label=false)
display(plt_OOP)
savefig(string(absPath,"/sweptPazyModalTipDispRangeDisp.pdf"))

# Tip twist vs. sweep angle for largest load
plt_twist = plot(xlabel="Sweep angle [deg]", ylabel="Tip twist [deg]")
plot!(rad2deg.(ΛRange), tip_twist[:,end], c=:black, m=:circle, lw=lw, ms=ms, label=false)
display(plt_twist)
savefig(string(absPath,"/sweptPazyModalTipDispRangeTwist.pdf"))

println("Finished sweptPazyModalTipDispRange.jl")