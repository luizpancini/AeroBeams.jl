using AeroBeams

# Wing pre-curvature range [1/m]
k2Range = [0.0, 0.018]

# Root pitch angle [rad]
θRange = deg2rad.(0:0.1:5)

# Airspeed range [m/s]
URange = vcat(0.5:0.5:40)

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
tip_u3 = Array{Float64}(undef,length(k2Range),length(θRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(k2Range),length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(k2Range),length(θRange))
flutterSpeed = Array{Float64}(undef,length(k2Range),length(θRange))
flutterFreq = Array{Float64}(undef,length(k2Range),length(θRange))
flutterTipDisp = Array{Float64}(undef,length(k2Range),length(θRange))

# Initialize guess solution
x0 = zeros(0)

# Sweep wing precurvature
for (i,k2) in enumerate(k2Range)
    # Sweep tip force
    for (j,θ) in enumerate(θRange)
        # Initialize guess solution
        x0 = zeros(0)
        # Sweep airspeed
        for (k,U) in enumerate(URange)
            # Display progress
            println("Solving for k2=$k2, θ = $(round(rad2deg(θ), digits=1)) deg, U = $U m/s")
            # Update model
            SMWFlutterPrecurvatureRange,_ = create_SMW(aeroSolver=aeroSolver,airfoil=airfoil,k2=k2,pitchAngle=θ,altitude=h,airspeed=U,nElem=nElem)
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
        freqs,damps,_ = mode_tracking(eigensols)
        # Get frequencies and dampings of oscillatory modes
        freqsOsc[i,j],dampsOsc[i,j],_ = filter_tracked_modes(freqs,damps; fmin=1e-2, maxTracks=nModes, tolDamp=1e-10)
        # Flutter speeds, frequencies and tip displacements
        flutterOnsetSpeeds,flutterOnsetFreqs,_,flutterOnsetTipDisp,_ = find_flutter(URange,freqsOsc[i,j],dampsOsc[i,j],tip_u3[i,j,:])
        # Set flutter speed and corresponding flutter frequency and tip displacement
        ind = isempty(flutterOnsetSpeeds) ? nothing : argmin(flutterOnsetSpeeds)
        flutterSpeed[i,j] = isnothing(ind) ? NaN : flutterOnsetSpeeds[ind]
        flutterFreq[i,j] = isnothing(ind) ? NaN : flutterOnsetFreqs[ind]
        flutterTipDisp[i,j] = isnothing(ind) ? NaN : flutterOnsetTipDisp[ind]
    end
end

## Plots
# Set paths
relPath = "/test/outputs/figures/SMWFlutterPrecurvatureRange2"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(k2Range), categorical=true)
lw = 2
ms = 3
gr()

# Initialize arrays
fSpeedOnsetMode = Array{Vector{Float64}}(undef,nModes,length(k2Range))
fSpeedOffsetMode = Array{Vector{Float64}}(undef,nModes,length(k2Range))

# Flutter speed vs. tip load
plt1 = plot(xlabel="Root angle [deg]", ylabel="Flutter speed [m/s]", xlims=[0,5], ylims=[0,40])
for (i,k2) in enumerate(k2Range)
    plot!([NaN], [NaN], c=colors[i], ls=:solid, lw=lw, label="k2=$k2 - onset")
    plot!([NaN], [NaN], c=colors[i], ls=:dash, lw=lw, label="k2=$k2 - offset")
end
for (i,k2) in enumerate(k2Range)
    plot!(rad2deg.(θRange), flutterSpeed[i,:], c=colors[i], ls=:solid, lw=lw, label=false)
end
display(plt1)
savefig(string(absPath,"/SMWFlutterPrecurvatureRange2_flutter.pdf"))

println("Finished SMWFlutterPrecurvatureRange2.jl")