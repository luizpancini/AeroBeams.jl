using AeroBeams

# Pitch angle range [deg]
θRange = deg2rad.([1.0, 2.2])

# Airspeed range [m/s]
URange = collect(0:1:40)

# Discretization
nElem = 20

# Number of modes to track
nModes = 4

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 50

# Initialize outputs
problem = Array{EigenProblem}(undef,length(θRange),length(URange))
eigensols = Array{EigenSolution}(undef,length(θRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(θRange))
tip_u3 = Array{Float64}(undef,length(θRange),length(URange))
tip_twist = Array{Float64}(undef,length(θRange),length(URange))

# Loop configurations
for (i,θ) in enumerate(θRange)
    # Model
    TDWingAirspeedRange = create_TDWing(pitchAngle=θ,nElem=nElem)
    # Initialize guess solution
    x0 = zeros(0)
    # Loop airspeeds
    for (j,U) in enumerate(URange)
        display("Solving for root angle = $(rad2deg(θ)) deg, U = $(round(Int,U)) m/s")
        # Update velocity of basis A 
        update_v_A!(TDWingAirspeedRange; v_A=[0,U,0])
        # Create and solve problem
        problem[i,j] = create_EigenProblem(model=TDWingAirspeedRange,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(problem[i,j])
        # Get outputs
        eigensols[i,j] = problem[i,j].eigensolution
        tipNodeFinalState = problem[i,j].loadstepHistory.states[end][nElem].nodal
        tip_u3[i,j] = tipNodeFinalState.u2[3]
        tip_twist[i,j] = twist_angle(tipNodeFinalState.p2_b)
        # Set initial guess for the next airspeed step
        x0 = problem[i,j].systemArrays.x
    end
    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols[i,:])
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[i],dampsOsc[i],_ = filter_tracked_modes(rads2Hz.(freqs),damps; fmin=eps(), fmax=40, maxTracks=nModes, dropOutsideBounds=true)
end

# Load reference solutions
u3_1deg_exp = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "u3_1deg_exp.txt"))
u3_1deg_num = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "u3_1deg_num.txt"))
u3_2_2deg_exp = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "u3_2_2deg_exp.txt"))
u3_2_2deg_num = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "u3_2_2deg_num.txt"))
th_1deg_exp = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "th_1deg_exp.txt"))
th_1deg_num = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "th_1deg_num.txt"))
th_2_2deg_exp = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "th_2_2deg_exp.txt"))
th_2_2deg_num = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "th_2_2deg_num.txt"))
freqs_ref = readdlm(pkgdir(AeroBeams,"test","referenceData","TDWing", "freqs.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/TDWingAirspeedRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Steady deformation at the last root angle and airspeed
deformationPlot = plot_steady_deformation(problem[end,end],view=(30,30),plotLegend=false,plotAxes=false,plotGrid=false,save=true,savePath=string(relPath,"/TDWingAirspeedRange_deformation.pdf"))
display(deformationPlot)

# Plot configurations
colors = cgrad(:darkrainbow, length(θRange), categorical=true)
colors2 = cgrad(:darkrainbow, 3, categorical=true)
lw = 2
ms = 4
msw = 0
labels = ["\\theta_{r} = $(rad2deg(θRange[1])) deg" "\\theta_{r} = $(rad2deg(θRange[2])) deg"]
gr()

# Tip flapwise displacement
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Tip flapwise displacement [m]")
plot!([NaN], [NaN], lc=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN], [NaN], lc=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
for i=eachindex(θRange)
    plot!(URange,tip_u3[i,:], c=colors[i], lw=lw, ls=:solid, label=false)
    if i==1
        plot!(u3_1deg_num[1,:],u3_1deg_num[2,:], c=colors[i], lw=lw, ls=:dash, label=false)
        scatter!(u3_1deg_exp[1,:],u3_1deg_exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        isfinite(tip_u3[i,2]) && annotate!(30, tip_u3[i,2], text(labels[i], :bottom, :left, colors[i]))
    else
        plot!(u3_2_2deg_num[1,:],u3_2_2deg_num[2,:], c=colors[i], lw=lw, ls=:dash, label=false)
        scatter!(u3_2_2deg_exp[1,:],u3_2_2deg_exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        isfinite(tip_u3[i,end-1]) && annotate!(35, tip_u3[i,end-1], text(labels[i], :top, :right, colors[i]))
    end
end
display(plt1)
savefig(string(absPath,"/TDWingAirspeedRange_disp.pdf"))

# Tip twist
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]")
plot!([NaN], [NaN], lc=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN], [NaN], lc=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!([NaN], [NaN ], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
for i=eachindex(θRange)
    plot!(URange,tip_twist[i,:], c=colors[i], lw=lw, ls=:solid, label=false)
    if i==1
        plot!(th_1deg_num[1,:],th_1deg_num[2,:], c=colors[i], lw=lw, ls=:dash, label=false)
        scatter!(th_1deg_exp[1,:],th_1deg_exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        isfinite(tip_twist[i,2]) && annotate!(30, tip_twist[i,2], text(labels[i], :bottom, :left, colors[i]))
    else
        plot!(th_2_2deg_num[1,:],th_2_2deg_num[2,:], c=colors[i], lw=lw, ls=:dash, label=false)
        scatter!(th_2_2deg_exp[1,:],th_2_2deg_exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        isfinite(tip_twist[i,end-1]) && annotate!(33, tip_twist[i,end-1], text(labels[i], :top, :right, colors[i]))
    end
end
display(plt2)
savefig(string(absPath,"/TDWingAirspeedRange_twist.pdf"))

# Aeroelastic frequencies at root angle of 1 deg
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]",legend=:outertop)
plot!([NaN], [NaN], lc=:black, lw=lw, ls=:solid, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="Arena & Lacarbonara (2013) - Num.")
modeLabels = ["Flapwise bending" "Chordwise bending" "Torsion"]
for (m,mode) in enumerate([1,2,4])
    plot!(URange, freqsOsc[1][mode,:], c=colors2[m], lw=lw, ls=:solid, label=false)
    scatter!(freqs_ref[1,:],freqs_ref[1+m,:], mc=colors2[m], ms=ms, msw=msw, label=false)
    isfinite(freqsOsc[1][mode,end]) && annotate!(10, freqsOsc[1][mode,end], text(modeLabels[m], :bottom, colors2[m]))
end
display(plt3)
savefig(string(absPath,"/TDWingAirspeedRange_freqs.pdf"))

println("Finished TDWingAirspeedRange.jl")