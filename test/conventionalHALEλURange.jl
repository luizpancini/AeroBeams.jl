using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Mode tracking option
modeTracking = true

# Aerodynamic solver
aeroSolver = Indicial()

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Discretization
nElemWing = 20
nElemTailBoom = 10
nElemHorzStabilizer = 10

# Set number of vibration modes
nModes = 15

# System solver for trim problem
relaxFactor = 0.5
maxIter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,displayStatus=false)

# Set stiffness and airspeed ranges, and initialize outputs
λRange = [1,5,50]
URange = collect(20:0.5:35)
untrackedFreqs = Array{Vector{Float64}}(undef,length(λRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(λRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(λRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(λRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(λRange),length(URange))
modeDampings = Array{Vector{Float64}}(undef,length(λRange),nModes)
modeFrequencies =  Array{Vector{Float64}}(undef,length(λRange),nModes)
eigenProblem = Array{EigenProblem}(undef,length(λRange),length(URange))

# Attachment springs
μ = 1e-1
ku = μ*[1; 1; 1]
kp = ku
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)
spring2 = create_Spring(elementsIDs=[nElemTailBoom],nodesSides=[2],ku=ku,kp=kp)

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving trim problem for λ = $λ, U = $U m/s")
        # Model for trim problem
        conventionalHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stiffnessFactor=λ,∞=1e12,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true)
        # Add springs
        add_springs_to_beam!(beam=tailBoom,springs=[spring1,spring2])
        # Update model
        conventionalHALEtrim.skipValidationMotionBasisA = true
        update_model!(conventionalHALEtrim)
        # Set initial guess solution as previous known solution
        x0Trim = j == 1 ? zeros(0) : trimProblem.x
        # Create and trim problem
        global trimProblem = create_TrimProblem(model=conventionalHALEtrim,systemSolver=NR,x0=x0Trim)
        solve!(trimProblem)
        # Extract trim variables
        trimAoA = trimProblem.aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
        trimThrust = stabilizersAero ? trimProblem.x[end-1]*trimProblem.model.forceScaling : trimProblem.x[end]*trimProblem.model.forceScaling
        trimδ = stabilizersAero ? trimProblem.x[end] : 0
        println("Trim AoA = $(trimAoA*180/π), trim thrust = $(trimThrust), trim δ = $(trimδ*180/π)")
        # Model for eigen problem
        conventionalHALEeigen,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stiffnessFactor=λ,∞=1e12,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=trimδ,thrust=trimThrust)
        # Add springs
        add_springs_to_beam!(beam=tailBoom,springs=[spring1,spring2])
        # Update model
        conventionalHALEeigen.skipValidationMotionBasisA = true
        update_model!(conventionalHALEeigen)
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=conventionalHALEeigen,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],jacobian=trimProblem.jacobian[1:end,1:end-trimProblem.model.nTrimVariables],inertia=trimProblem.inertia)
        solve_eigen!(eigenProblem[i,j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = eigenProblem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(eigenProblem[i,j].dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = eigenProblem[i,j].eigenvectorsOscillatoryCplx
    end
    # Frequencies and dampings after mode tracking
    if modeTracking
        freqs[i,:],damps[i,:],_ = mode_tracking(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    else
        freqs[i,:],damps[i,:] = untrackedFreqs[i,:],untrackedDamps[i,:]
    end
    # Separate frequencies and dampings by mode
    for mode in 1:nModes
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URange)]
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URange)]
    end
end

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(λRange)))
lw = 2
ms = 5
lstyle = [:solid, :dash, :dot]
mshape = [:circle, :star, :utriangle]
relPath = "/test/outputs/figures/conventionalHALEλURange"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Mode shapes of flexible aircraft at lowest airspeed
modesPlot = plot_mode_shapes(eigenProblem[1,1],nModes=5,scale=10,view=(30,30),legendPos=:outertop,save=true,savePath=string(relPath,"/conventionalHALEλURange_modeShapes.pdf"))
display(modesPlot)
# Root locus
gr()
plt0 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", ylims=[0,50])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=0, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=0, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=0, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
display(plt0)
savefig(string(absPath,"/conventionalHALEλURange_rootlocus.pdf"))
# Root locus
plt1 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,0], ylims=[0,10])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=0, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=0, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=0, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/conventionalHALEλURange_rootlocus2.pdf"))
# Root locus (zoom)
plt2 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.1,0.1], ylims=[0,0.5])
scatter!([NaN], [NaN], c=colors[1], shape=mshape[1], ms=ms, msw=0, label="Flexible")
scatter!([NaN], [NaN], c=colors[2], shape=mshape[2], ms=ms, msw=0, label="Stiffened")
scatter!([NaN], [NaN], c=colors[3], shape=mshape[3], ms=ms, msw=0, label="Rigid")
for (i,λ) in enumerate(λRange)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
display(plt2)
savefig(string(absPath,"/conventionalHALEλURange_rootlocuszoom.pdf"))

println("Finished conventionalHALEλURange.jl")