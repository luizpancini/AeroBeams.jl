using AeroBeams, LinearAlgebra, LinearInterpolations, Plots, ColorSchemes, DelimitedFiles

# Option for mode tracking
modeTracking = true

# Pazy wing
wing,L,nElem,chord,normSparPos,airfoil,surf = create_Pazy(p0=[0;-π/2;0], airfoil=flatPlate)

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
PazyWingFlutterTipMassRange = create_Model(name="PazyWingFlutterTipMassRange",beams=[wing],BCs=[clamp],gravityVector=[0;0;-9.807])

# Set system solver options (limit initial load factor)
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Number of modes
nModes = 3

# Set tip mass and its position ranges, airspeed range, and initialize outputs
configurations = [1; 2; 3]
tipMassRange = [0; 0.01; 0.01]
tipMassPosRange = chord/2*[0; 1; -1]
URange = collect(0:0.5:100)
untrackedFreqs = Array{Vector{Float64}}(undef,3,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,3,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,3,length(URange))
freqs = Array{Vector{Float64}}(undef,3,length(URange))
damps = Array{Vector{Float64}}(undef,3,length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,3,nModes)
modeDampings = Array{Vector{Float64}}(undef,3,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,3,nModes)

# Sweep configurations
for c in configurations
    # Reset point inertias on beam
    wing.pointInertias = Vector{PointInertia}()
    update_beam!(wing)
    # Create tip mass 
    tipMass = PointInertia(elementID=nElem,η=[L/nElem/2;tipMassPosRange[c];0],mass=tipMassRange[c])
    # Add tip mass to beam
    add_point_inertias_to_beam!(wing,inertias=[tipMass])
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for configuration $c, U = $U m/s")
        # Set tip loss function at current airspeed and root angle
        surf.tipLossDecayFactor = tip_loss_factor_Pazy(0,U)
        update_beam!(wing)
        # Update velocity of basis A (and update model)
        set_motion_basis_A!(model=PazyWingFlutterTipMassRange,v_A=[0;U;0])
        # Create and solve problem
        problem = create_EigenProblem(model=PazyWingFlutterTipMassRange,nModes=nModes,systemSolver=NR)
        solve!(problem)
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[c,j] = problem.frequenciesOscillatory
        untrackedDamps[c,j] = round_off!(problem.dampingsOscillatory,1e-8)
        untrackedEigenvectors[c,j] = problem.eigenvectorsOscillatoryCplx
    end
    # Apply mode tracking, if applicable
    if modeTracking
        freqs[c,:],damps[c,:],_ = mode_tracking(URange,untrackedFreqs[c,:],untrackedDamps[c,:],untrackedEigenvectors[c,:])
    else
        freqs[c,:],damps[c,:] = untrackedFreqs[c,:],untrackedDamps[c,:]
    end
    # Separate frequencies and damping ratios by mode
    for mode in 1:nModes
        modeFrequencies[c,mode] = [freqs[c,i][mode] for i in eachindex(URange)]
        modeDampings[c,mode] = [damps[c,i][mode] for i in eachindex(URange)]
        modeDampingRatios[c,mode] = modeDampings[c,mode]./modeFrequencies[c,mode]
    end
end


# Plots
# ------------------------------------------------------------------------------
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
lstyles = [:solid :dash :dot]
relPath = "/test/outputs/figures/PazyWingFlutterTipMassRange"
absPath = string(pwd(),relPath)
mkpath(absPath)
gr()
# V-g-f for all configurations
plt11 = plot(ylabel="Frequency [Hz]")
for c in configurations
    for mode in 1:nModes
        plot!(URange, modeFrequencies[c,mode]/(2*π), c=modeColors[mode], ls=lstyles[c], lw=lw, label=false)
    end
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.1,0.05], legend=:topleft)
plot!([NaN], [NaN], c=:black, ls=lstyles[1], lw=lw, label="Clean wing")
plot!([NaN], [NaN], c=:black, ls=lstyles[2], lw=lw, label="LE weight")
plot!([NaN], [NaN], c=:black, ls=lstyles[3], lw=lw, label="TE weight")
for c in configurations
    for mode in 1:nModes
        plot!(URange, modeDampingRatios[c,mode], c=modeColors[mode], ls=lstyles[c], lw=lw, label=false)
    end
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/PazyWingFlutterTipMassRange.pdf"))

println("Finished PazyWingFlutterTipMassRange.jl")