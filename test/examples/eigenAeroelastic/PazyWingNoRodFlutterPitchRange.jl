using AeroBeams

# Pitch angle range [rad]
θRange = unique([vcat(0:0.5:1)...,vcat(1:1:7)...])

# Airspeed range [m/s]
URange = vcat(1e-3,1:1:130)

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(FlatPlate)

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "Exponential"

# Flag for upright position
upright = true

# Mesh: :nominal or :refined
mesh = :refined

# Fixed geometrical and discretization properties
nElem,L,chord,normSparPos,_ = geometrical_properties_Pazy(; mesh=mesh)

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = mesh == :nominal ? 45 : 65

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(θRange),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(θRange))
dampsOsc = Array{Matrix{Float64}}(undef,length(θRange))
flutterOnsetSpeeds = Array{Vector{Float64}}(undef,length(θRange))
flutterOnsetFreqs = Array{Vector{Float64}}(undef,length(θRange))
flutterOnsetModes = Array{Vector{Int}}(undef,length(θRange))
flutterOffsetSpeeds = Array{Vector{Float64}}(undef,length(θRange))
flutterOffsetFreqs = Array{Vector{Float64}}(undef,length(θRange))
flutterOffsetModes = Array{Vector{Int}}(undef,length(θRange))

# Sweep root angle
for (i,θ) in enumerate(θRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(θ, digits=2)) deg, U = $U m/s")
        # Model
        PazyWingNoRodFlutterPitchRange,_ = create_Pazy(aeroSolver=aeroSolver,withWingtipRod=false,airfoil=airfoil,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,upright=upright,pitchAngle=deg2rad(θ),airspeed=U,mesh=mesh)
        # Create and solve problem
        problem = create_EigenProblem(model=PazyWingNoRodFlutterPitchRange,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(problem)
        # Eigensolution
        eigensols[i,j] = problem.eigensolution
        # Set initial guess solution for next iteration
        x0 = problem.systemArrays.x
    end
    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols[i,:])
    rads2Hz!(freqs)
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[i],dampsOsc[i],_ = filter_tracked_modes(freqs,damps; fmin=1e0, maxTracks=nModes)
    # Flutter onset and offset speeds, frequencies and tip OOP displacements
    flutterOnsetSpeeds[i],flutterOnsetFreqs[i],flutterOnsetModes[i],flutterOffsetSpeeds[i],flutterOffsetFreqs[i],flutterOffsetModes[i] = find_flutter(URange,freqsOsc[i],dampsOsc[i])
end

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingNoRodFlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ts = 10
fs = 16
lfs = 9
lw = 2
ms = 10
ms2 = 4
msw = 0
gr()

# Extract flutter onset and offset data for the hump mode
humpMode = 3
UfOn = Float64[]
UfOff = Float64[]
for i in eachindex(θRange)
    for m in eachindex(flutterOnsetModes[i])
        if flutterOnsetModes[i][m] == humpMode
            push!(UfOn, flutterOnsetSpeeds[i][m])
            break
        end
    end
    for m in eachindex(flutterOffsetModes[i])
        if flutterOffsetModes[i][m] == humpMode
            push!(UfOff, flutterOffsetSpeeds[i][m])
            break
        end
    end
end

# Flutter onset and offset speeds vs root pitch angle
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]", xlims=[30,125], ylims=[0,7.25], xticks=collect(30:10:125), yticks=collect(0:1:7), legend=:topright, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!(Shape(vcat(UfOn,reverse(UfOff)),vcat(θRange,reverse(θRange))), fillcolor = plot_color(:red, 0.25), lw=lw, label="Hump flutter region")
display(plt1)
savefig(string(absPath,"/PazyWingNoRodFlutterPitchRange_flutterBoundaryPitch.pdf"))

println("Finished PazyWingNoRodFlutterPitchRange.jl")