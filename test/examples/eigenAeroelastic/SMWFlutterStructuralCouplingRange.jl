using AeroBeams

# Torsion / IP bending coupling factor range
ΨRange = [-0.2, 0, 0.2]

# Tip load range
F3Range = collect(0:1:35)

# Airspeed range
URange = collect(20:0.5:35)

# Aerodynamic solver
aeroSolver = Indicial()

# Altitude 
h = 20e3

# Gravity
g = 0

# Discretization
nElem = 16

# Number of modes to track
nModes = 10

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 50

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(ΨRange),length(F3Range),length(URange))
flutterOnsetSpeeds = Array{Vector{Float64}}(undef,length(ΨRange),length(F3Range))
flutterSpeed = Array{Float64}(undef,length(ΨRange),length(F3Range))

# Sweep structural coupling 
for (i,Ψ) in enumerate(ΨRange)
    # Sweep tip force
    for (j,F3) in enumerate(F3Range)
        # Initialize guess solution
        x0 = zeros(0)
        # Sweep airspeed
        for (k,U) in enumerate(URange)
            # Display progress
            println("Solving for Ψ=$Ψ, F=$F3 N, U = $U m/s")
            # Create model
            model,_ = create_SMW(aeroSolver=aeroSolver,nElem=nElem,altitude=h,airspeed=U,g=g,Ψ=Ψ,tipF3=F3)
            # Create and solve problem
            problem = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev,x0=x0)
            solve!(problem)
            # Eigensolution
            eigensols[i,j,k] = problem.eigensolution
            # Set initial guess solution for next iteration
            x0 = problem.systemArrays.x
        end
        # Apply mode tracking
        freqs,damps,_ = mode_tracking(eigensols[i,j,:])
        # Get frequencies and dampings of oscillatory modes
        freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs,damps; fmin=1, maxTracks=nModes, tolDamp=1e-9)
        # Flutter speed
        flutterOnsetSpeeds[i,j],_ = find_flutter(URange,freqsOsc,dampsOsc)
        flutterSpeed[i,j] = flutterOnsetSpeeds[i,j][1]
    end
end

# Load reference data
flutterSpeedVsTipLoadΨ0 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "flutterSpeedVsTipLoadPsi0_0.txt"))
flutterSpeedVsTipLoadΨp02 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "flutterSpeedVsTipLoadPsi0_2.txt"))
flutterSpeedVsTipLoadΨm02 = readdlm(pkgdir(AeroBeams,"test","referenceData","SMW", "flutterSpeedVsTipLoadPsi-0_2.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/SMWFlutterStructuralCouplingRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(ΨRange), categorical=true)
lw = 2
ms = 3
gr()

# Flutter speed vs tip load for several structural couplings
plt1 = plot(xlabel="Tip Load [N]", ylabel="Flutter speed [m/s]", xlims=[0,35], ylims=[20,35])
plot!([NaN], [NaN], c=:black, lw=lw, legend=:bottomleft, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Patil et al. (2001)")
for (i,Ψ) in enumerate(ΨRange)
    plot!(F3Range, flutterSpeed[i,:], c=colors[i], lw = lw, label="\\Psi=$Ψ")
    if Ψ == -0.2
        scatter!(flutterSpeedVsTipLoadΨm02[1,:], flutterSpeedVsTipLoadΨm02[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif Ψ == 0.0
        scatter!(flutterSpeedVsTipLoadΨ0[1,:], flutterSpeedVsTipLoadΨ0[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif Ψ == +0.2
        scatter!(flutterSpeedVsTipLoadΨp02[1,:], flutterSpeedVsTipLoadΨp02[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/SMWFlutterStructuralCouplingRange.pdf"))

println("Finished SMWFlutterStructuralCouplingRange.jl")