using AeroBeams

# Root angle range [rad]
θRange = deg2rad.(0:0.2:5.0)

# Airspeed range [m/s]
URange = vcat(1e-3,0.25:0.25:35)

# Altitude [m]
h = 20e3

# Discretization
nElem = 16

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil
airfoil = deepcopy(FlatPlate)

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 60

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(URange))
tip_u3 = Array{Float64}(undef,length(θRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(θRange))
flutterOnsetSpeeds = Array{Vector{Float64}}(undef,length(θRange))
flutterOnsetFreqs = Array{Vector{Float64}}(undef,length(θRange))
flutterOnsetModes = Array{Vector{Int}}(undef,length(θRange))
flutterOnsetTipDisp = Array{Vector{Float64}}(undef,length(θRange))
flutterOffsetSpeeds = Array{Vector{Float64}}(undef,length(θRange))
flutterOffsetFreqs = Array{Vector{Float64}}(undef,length(θRange))
flutterOffsetModes = Array{Vector{Int}}(undef,length(θRange))
flutterOffsetTipDisp = Array{Vector{Float64}}(undef,length(θRange))

# Initialize guess solution
x0 = zeros(0)

# Sweep root angle
for (i,θ) in enumerate(θRange)
    # Update root angle
    SMWFlutterPitchRange,_ = create_SMW(aeroSolver=aeroSolver,airfoil=airfoil,pitchAngle=θ,altitude=h,nElem=nElem)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(rad2deg(θ), digits=1)) deg, U = $U m/s")
        # Update airspeed
        update_v_A!(SMWFlutterPitchRange,v_A=[0,U,0])
        # Create and solve problem
        problem = create_EigenProblem(model=SMWFlutterPitchRange,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(problem)
        # Eigensolution
        eigensols[j] = problem.eigensolution
        # Tip OOP displacement
        tip_u3[i,j] = problem.loadstepHistory.states[end][nElem].nodal.u2_b[3]
        # Set initial guess solution for next iteration
        x0 = problem.systemArrays.x
    end
    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols)
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[i],dampsOsc[i],_ = filter_tracked_modes(freqs,damps; fmin=1e-2, maxTracks=nModes)
    # Flutter speeds, frequencies and tip displacements at current root angle
    flutterOnsetSpeeds[i],flutterOnsetFreqs[i],flutterOnsetModes[i],flutterOnsetTipDisp[i],flutterOffsetSpeeds[i],flutterOffsetFreqs[i],flutterOffsetModes[i],flutterOffsetTipDisp[i] = find_flutter(URange,freqsOsc[i],dampsOsc[i],tip_u3[i,:])
end

# Load reference data
flutterSpeedRef = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "flutterSpeedVsRootAoA.txt"))
flutterFreqRef = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "flutterFreqVsRootAoA.txt"))
flutterTipDispRef = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "flutterTipDispVsRootAoA.txt"))
speedVsDispRootAoA0 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "speedVsDispRootAoA0_0.txt"))
speedVsDispRootAoA05 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "speedVsDispRootAoA0_5.txt"))
speedVsDispRootAoA1 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "speedVsDispRootAoA1_0.txt"))
speedVsDispRootAoA2 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "speedVsDispRootAoA2_0.txt"))
freqVsSpeedRootAoA2 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "V-f-RootAoA2.txt"))
dampVsSpeedRootAoA2 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "V-g-RootAoA2.txt"))
freqVsSpeedRootAoA3 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "V-f-RootAoA3.txt"))
dampVsSpeedRootAoA3 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "V-g-RootAoA3.txt"))
freqVsSpeedRootAoA5 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "V-f-RootAoA5.txt"))
dampVsSpeedRootAoA5 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "V-g-RootAoA5.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/SMWFlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
markers = [:utriangle, :dtriangle]
lw = 2
ms = 3
msw = 0
flutterOnsetModesLabels = ["1st OOP - onset", "2nd OOP - onset", "1st T-IP - onset", "3rd OOP - onset", "2nd T-IP - onset"]
flutterOffsetModesLabels = ["1st OOP - offset", "2nd OOP - offset", "1st T-IP - offset", "3rd OOP - offset", "2nd T-IP - offset"]
gr()

# Helper function to plot flutter data over root pitch angle
flutterModes = sort(unique(reduce(vcat, flutterOnsetModes)))
function plot_overθ(plt,flutterModes,θRange,modes,param,labels,marker)
    for m in flutterModes
        xs = Float64[]
        ys = Float64[]
        for i in eachindex(θRange)
            for (mj, Uj) in zip(modes[i], param[i])
                if mj == m
                    push!(xs, rad2deg(θRange[i]))
                    push!(ys, Uj)
                end
            end
        end
        scatter!(plt, xs, ys, marker=marker, c=modeColors[m], ms=ms, msw=msw, label=labels[m])
    end
end

# Flutter speed and frequency vs. root pitch angle
plt11 = plot(ylabel="Flutter speed [m/s]", xlims=[0,5], ylims=[0,35], legend=:bottomleft)
plot_overθ(plt11, flutterModes, θRange, flutterOnsetModes, flutterOnsetSpeeds, flutterOnsetModesLabels, markers[1])
plot_overθ(plt11, flutterModes, θRange, flutterOffsetModes, flutterOffsetSpeeds, flutterOffsetModesLabels, markers[2])
plot!(flutterSpeedRef[1,:], flutterSpeedRef[2,:], c=:black, ls=:solid, lw=lw, label="Patil et al. (2001)")

plt12 = plot(xlabel="Root pitch angle, \$\\theta\$ [deg]", ylabel="Flutter frequency [rad/s]", xlims=[0,5], ylims=[0,35])
plot_overθ(plt12, flutterModes, θRange, flutterOnsetModes, flutterOnsetFreqs, flutterOnsetModesLabels, markers[1])
plot_overθ(plt12, flutterModes, θRange, flutterOffsetModes, flutterOffsetFreqs, flutterOffsetModesLabels, markers[2])
plot!(flutterFreqRef[1,:], flutterFreqRef[2,:], c=:black, ls=:solid, lw=lw, label=false)

plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/SMWFlutterPitchRange_flutterSpeedAndFreq.pdf"))

# Flutter tip displacement vs. root pitch angle
plt2 = plot(xlabel="Root pitch angle, \$\\theta\$ [deg]",ylabel="Flutter tip displacement [m]", xlims=[0,5], ylims=[-3,3])
plot_overθ(plt2, flutterModes, θRange, flutterOnsetModes, flutterOnsetTipDisp, flutterOnsetModesLabels, markers[1])
plot_overθ(plt2, flutterModes, θRange, flutterOffsetModes, flutterOffsetTipDisp, flutterOffsetModesLabels, markers[2])
plot!(flutterTipDispRef[1,:], flutterTipDispRef[2,:], c=:black, ls=:solid, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/SMWFlutterPitchRange_flutterDisp.pdf"))

# Tip displacement vs airspeed
θplot = deg2rad.([0, 0.5, 1.0, 2.0])
colors = cgrad(:rainbow, length(θplot), categorical=true)
indθ2plot = findall(vec(any(θRange .== θplot', dims=2)))
plt3 = plot(xlabel="Tip displacement [m]",ylabel="Airspeed [m/s]", xlims=[-3,3], ylims=[0,35])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="Patil et al. (2001)")
for i in eachindex(indθ2plot)
    plot!(tip_u3[indθ2plot[i],:], URange, c=colors[i], lw=lw, label="\\theta=$(rad2deg(θplot[i])) deg")
    if i == 1
        scatter!(speedVsDispRootAoA0[1,:],speedVsDispRootAoA0[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif i == 2
        scatter!(speedVsDispRootAoA05[1,:],speedVsDispRootAoA05[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif i == 3
        scatter!(speedVsDispRootAoA1[1,:],speedVsDispRootAoA1[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif i == 4
        scatter!(speedVsDispRootAoA2[1,:],speedVsDispRootAoA2[2,:], mc=colors[i], ms=ms, msw=0, label=false)    
    end
end
display(plt3)
savefig(string(absPath,"/SMWFlutterPitchRange_tipDisp.pdf"))

# V-g-f of selected mode for selected root angles
mode2plot = 2
θplot = deg2rad.([2.0, 3.0, 5.0])
colors = cgrad(:rainbow, length(θplot), categorical=true)
indθ2plot = findall(vec(any(θRange .== θplot', dims=2)))
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[10,30], ylims=[0,35], legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="Patil et al. (2001)")
for i in eachindex(indθ2plot)
    plot!(URange, freqsOsc[indθ2plot[i]][mode2plot,:], c=colors[i], lw=lw, label=false)
    if i == 1
        scatter!(freqVsSpeedRootAoA2[1,:],freqVsSpeedRootAoA2[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif i == 2
        scatter!(freqVsSpeedRootAoA3[1,:],freqVsSpeedRootAoA3[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif i == 3
        scatter!(freqVsSpeedRootAoA5[1,:],freqVsSpeedRootAoA5[2,:], mc=colors[i], ms=ms, msw=0, label=false)   
    end
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", xlims=[10,30], ylims=[-1.,0.25], legend=:bottomleft)
for i in eachindex(indθ2plot)
    plot!(URange, dampsOsc[indθ2plot[i]][mode2plot,:], c=colors[i], lw=lw, label="\$\\theta=$(round(Int, rad2deg(θplot[i])))\$ deg")
    if i == 1
        scatter!(dampVsSpeedRootAoA2[1,:],dampVsSpeedRootAoA2[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif i == 2
        scatter!(dampVsSpeedRootAoA3[1,:],dampVsSpeedRootAoA3[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif i == 3
        scatter!(dampVsSpeedRootAoA5[1,:],dampVsSpeedRootAoA5[2,:], mc=colors[i], ms=ms, msw=0, label=false)   
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,"/SMWFlutterPitchRange_Vgf.pdf"))

println("Finished SMWFlutterPitchRange.jl")