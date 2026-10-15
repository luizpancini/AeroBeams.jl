using AeroBeams

# Wing pre-curvature range [1/m]
k2Range = [0.0, 0.018]

# Tip force range [N]
F3Range = collect(0:2:40)

# Airspeed range [m/s]
URange = collect(20:0.5:35)

# Altitude [m]
h = 20e3

# Pitch angle [rad]
θ = deg2rad(0)

# Gravity
g = 0

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
nev = 50

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(URange))
tip_u3 = Array{Float64}(undef,length(k2Range),length(F3Range),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(k2Range),length(F3Range))
dampsOsc = Array{Matrix{Float64}}(undef,length(k2Range),length(F3Range))
flutterSpeed = Array{Float64}(undef,length(k2Range),length(F3Range))
flutterFreq = Array{Float64}(undef,length(k2Range),length(F3Range))
flutterTipDisp = Array{Float64}(undef,length(k2Range),length(F3Range))

# Initialize guess solution
x0 = zeros(0)

# Sweep wing precurvature
for (i,k2) in enumerate(k2Range)
    # Sweep tip force
    for (j,F3) in enumerate(F3Range)
        # Initialize guess solution
        x0 = zeros(0)
        # Sweep airspeed
        for (k,U) in enumerate(URange)
            # Display progress
            println("Solving for k2=$k2, F3 = $F3 N, U = $U m/s")
            # Update model
            SMWFlutterPrecurvatureRange,_ = create_SMW(aeroSolver=aeroSolver,airfoil=airfoil,k2=k2,pitchAngle=θ,altitude=h,g=g,tipF3=F3,airspeed=U,nElem=nElem)
            # Create and solve problem
            problem = create_EigenProblem(model=SMWFlutterPrecurvatureRange,eigenmethod=eigenmethod,nev=nev,x0=x0)
            solve!(problem)
            # Eigensolution
            eigensols[k] = problem.eigensolution
            # Tip OOP displacement
            tip_u3[i,j,k] = problem.loadstepHistory.states[end][nElem].nodal.u2_b[3]
            # Set initial guess solution for next iteration
            x0 = problem.systemArrays.x
        end
        # Apply mode tracking
        freqs,damps,_ = mode_tracking(eigensols; α=1e0)
        # Get frequencies and dampings of oscillatory modes
        freqsOsc[i,j],dampsOsc[i,j],_ = filter_tracked_modes(freqs,damps; fmin=1e-2, maxTracks=nModes, tolDamp=1e-9)
        # Flutter speeds, frequencies and tip displacements
        flutterOnsetSpeeds,flutterOnsetFreqs,_,flutterOnsetTipDisp,_ = find_flutter(URange,freqsOsc[i,j],dampsOsc[i,j],tip_u3[i,j,:])
        # Set flutter speed and corresponding flutter frequency and tip displacement
        flutterSpeed[i,j] = flutterOnsetSpeeds[1]
        flutterFreq[i,j] = flutterOnsetFreqs[1]
        flutterTipDisp[i,j] = flutterOnsetTipDisp[1]
    end
end

# Load reference data
flutterSpeedVsTipLoadk0 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "flutterSpeedVsTipLoadk0.txt"))
flutterSpeedVsTipLoadk2 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "flutterSpeedVsTipLoadk2.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/SMWFlutterPrecurvatureRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(k2Range), categorical=true)
modeColors = cgrad(:rainbow, nModes, categorical=true)
mshape = [:utriangle, :utriangle, :circle, :circle, :utriangle]
lw = 2
ms = 4
gr()

# Flutter speed vs. tip load
plt1 = plot(xlabel="Tip load [N]", ylabel="Flutter speed [m/s]", xlims=[0,40], ylims=[0,40])
for (i,k2) in enumerate(k2Range)
    plot!(F3Range, flutterSpeed[i,:], c=colors[i], lw=lw, label="AeroBeams - \$k_2\$=$k2")
end
plot!(flutterSpeedVsTipLoadk0[1,:], flutterSpeedVsTipLoadk0[2,:], c=:black, ls=:dash, lw=lw, label="Patil et al. (2001)")
plot!(flutterSpeedVsTipLoadk2[1,:], flutterSpeedVsTipLoadk2[2,:], c=:black, ls=:dash, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/SMWFlutterPrecurvatureRange_flutter.pdf"))

# Flutter speed vs. tip disp
plt2 = plot(xlabel="Tip OOP displacement [% semispan]", ylabel="Flutter speed [m/s]", xlims=[0,20], ylims=[0,40])
for (i,k2) in enumerate(k2Range)
    plot!(flutterTipDisp[i,:]/16*100, flutterSpeed[i,:], c=colors[i], lw=lw, label="\$k_2\$=$k2")
end
display(plt2)
savefig(string(absPath,"/SMWFlutterPrecurvatureRange_flutterTipDisp.pdf"))

# Frequency vs tip displacement of initially straight wing (at highest load)
Uind = 1
plt3 = plot(xlabel = "Tip OOP displacement [m]", ylabel="Frequency [rad/s]", ylims=[0,50])
scatter!([NaN], [NaN], c=:black, shape=mshape[1], ms=ms, msw=msw, label="OOP bending modes")
scatter!([NaN], [NaN], c=:black, shape=mshape[3], ms=ms, msw=msw, label="T-IP bending modes")
for j in eachindex(F3Range)
    for mode in 1:nModes
        scatter!([tip_u3[1,j,Uind]], [freqsOsc[1,j][mode,Uind]], c=modeColors[mode], shape=mshape[mode], ms=ms, msw=msw, label=false)
    end
end
display(plt3)
savefig(string(absPath,"/SMWFlutterPrecurvatureRange_freq_vs_disp.pdf"))

println("Finished SMWFlutterPrecurvatureRange.jl")