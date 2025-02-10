using AeroBeams

# Wing airfoil
wingAirfoil = deepcopy(NACA23012A)

# Option for reduced chord
reducedChord = true

# Option to include beam pods
beamPods = true

# Option to set payload on wing
payloadOnWing = true

# Aerodynamic solvers
aeroSolvers = [Indicial(), BLi()]

# Payload range
PRange = collect(0:10:500)

# Airspeed [m/s]
U = 40*0.3048

# Number of modes
nModes = 10

# System solver for trim problem
relaxFactor = 0.5
maxIter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter)

# Attachment springs
μ = 1e-2
ku = μ*[1; 1; 1]
kp = ku
spring = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(aeroSolvers),length(PRange))
eigenProblem = Array{EigenProblem}(undef,length(aeroSolvers),length(PRange))
trimAoA = Array{Float64}(undef,length(aeroSolvers),length(PRange))
trimThrust = Array{Float64}(undef,length(aeroSolvers),length(PRange))
trimδ = Array{Float64}(undef,length(aeroSolvers),length(PRange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(aeroSolvers),length(PRange))
untrackedDamps = Array{Vector{Float64}}(undef,length(aeroSolvers),length(PRange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(aeroSolvers),length(PRange))
freqs = Array{Vector{Float64}}(undef,length(aeroSolvers),length(PRange))
damps = Array{Vector{Float64}}(undef,length(aeroSolvers),length(PRange))
modeFrequencies = Array{Vector{Float64}}(undef,length(aeroSolvers),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(aeroSolvers),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(aeroSolvers),nModes)

# Sweep aero solvers
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Sweep payload
    for (j,P) in enumerate(PRange)
        # Display progress
        println("Solving for aeroSolver $i, payload = $P lb")
        # Model for trim problem
        heliosFlutter,midSpanElem,_,_,rightWingStraight,_ = create_Helios(reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,aeroSolver=aeroSolver,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
        # Add springs at wing root
        add_springs_to_beam!(beam=rightWingStraight,springs=[spring])
        # Update model
        heliosFlutter.skipValidationMotionBasisA = true
        update_model!(heliosFlutter)
        # Set initial guess solution as previous known solution
        x0Trim = (j==1) ? zeros(0) : trimProblem[i,j-1].x
        # Create and solve trim problem
        trimProblem[i,j] = create_TrimProblem(model=heliosFlutter,systemSolver=NR,x0=x0Trim)
        solve!(trimProblem[i,j])
        # Extract trim variables
        trimAoA[i,j] = trimProblem[i,j].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.αₑ*180/π
        trimThrust[i,j] = trimProblem[i,j].x[end-1]*trimProblem[i,j].model.forceScaling
        trimδ[i,j] = trimProblem[i,j].x[end]
        println("Trim summary: AoA = $(trimAoA[i,j]), T = $(trimThrust[i,j]), δ = $(trimδ[i,j]*180/π)")
        # Model for eigen problem
        heliosEigen,_,_,_,rightWingStraight,_ = create_Helios(reducedChord=reducedChord,payloadOnWing=payloadOnWing,beamPods=beamPods,wingAirfoil=wingAirfoil,aeroSolver=aeroSolver,payloadPounds=P,airspeed=U,δ=trimδ[i,j],thrust=trimThrust[i,j])
        # Add springs at wing root
        add_springs_to_beam!(beam=rightWingStraight,springs=[spring])
        # Update model
        heliosEigen.skipValidationMotionBasisA = true
        update_model!(heliosEigen)
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=heliosEigen,nModes=nModes,frequencyFilterLimits=[5e-3,Inf64],jacobian=trimProblem[i,j].jacobian[1:end,1:end-2],inertia=trimProblem[i,j].inertia,refTrimProblem=trimProblem[i,j])
        solve_eigen!(eigenProblem[i,j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = eigenProblem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(eigenProblem[i,j].dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = eigenProblem[i,j].eigenvectorsOscillatoryCplx
    end
    # Apply mode tracking, if applicable
    freqs[i,:],damps[i,:],_ = mode_tracking(PRange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies and damping ratios by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,k][mode] for k in eachindex(PRange)]
        modeDampings[i,mode] = [damps[i,k][mode] for k in eachindex(PRange)]
        modeDampingRatios[i,mode] = modeDampings[i,mode]./modeFrequencies[i,mode]
    end
end

# Set paths
relPath = "/dev/helios/figures/heliosFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
using Plots, ColorSchemes
gr()
ts = 10
fs = 16
lw = 2
ms = 3
msw = 0
solversNames = ["Indical", "BLi"]
solversLabels = ["Linear", "Dynamic stall"]
solversLS = [:dash, :solid]
modeColors = cgrad(:rainbow, nModes, categorical=true)

# Trim root angle of attack
plt_trimAoA = plot(xlabel="Payload [lb]", ylabel="Trim root AoA [deg]", xlims=[0,500], ylims=[0,6], tickfont=font(ts), guidefont=font(fs), legend=:bottomright, legendfontsize=12)
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(PRange, trimAoA[i,:], c=:black, lw=lw, ls=solversLS[i], label=solversLabels[i])
end
display(plt_trimAoA)
savefig(string(absPath,"/heliosFlutter_trimAoA.pdf"))

# Trim elevator deflection
plt_trimDelta = plot(xlabel="Payload [lb]", ylabel="Trim elevator deflection [deg]", xlims=[0,500], ylims=[0,12], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(PRange, trimδ[i,:]*180/π, c=:black, lw=lw, ls=solversLS[i], label=false)
end
display(plt_trimDelta)
savefig(string(absPath,"/heliosFlutter_trimDelta.pdf")) 

# Trim thrust per motor
plt_trimThrust = plot(xlabel="Payload [lb]", ylabel="Trim thrust per motor [N]", xlims=[0,500], ylims=[0,Inf], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(PRange, trimThrust[i,:], c=:black, lw=lw, ls=solversLS[i], label=false)
end
display(plt_trimThrust)
savefig(string(absPath,"/heliosFlutter_trimThrust.pdf")) 

# V-g-f
for (i,aeroSolver) in enumerate(aeroSolvers)
    plt_Vf = plot(ylabel="Frequency [rad/s]", tickfont=font(ts), guidefont=font(14), ylims=[0,10])
    for mode in 1:nModes
        scatter!(PRange, modeFrequencies[i,mode], c=modeColors[mode], ms=ms, msw=0, label=false)
    end
    plt_Vg = plot(xlabel="Payload [lb]", ylabel="Damping ratio", ylims=[-1,0.25], tickfont=font(ts), guidefont=font(14))
    for mode in 1:nModes
        scatter!(PRange, modeDampingRatios[i,mode], c=modeColors[mode], ms=ms, msw=msw, label=false)
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(string(absPath,"/heliosFlutter_Vgf_",solversNames[i],".pdf"))
end

# Root locus
for (i,aeroSolver) in enumerate(aeroSolvers)
    plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", tickfont=font(ts), guidefont=font(fs), xlims=[-5,0.5], ylims=[0,10])
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=modeColors[mode], ms=ms, msw=msw, label=false)
    end
    display(plt_RL)
    savefig(string(absPath,"/heliosFlutter_rootlocus_",solversNames[i],".pdf"))
end

# Root locus (phugoid zoom)
for (i,aeroSolver) in enumerate(aeroSolvers)
    plt_phugoid = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.1,0.2], ylims=[0,0.6], tickfont=font(ts), guidefont=font(fs))
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=modeColors[mode], ms=ms, msw=msw, label=false)
    end
    display(plt_phugoid)
    savefig(string(absPath,"/heliosFlutter_Phzoom_",solversNames[i],".pdf"))
end

println("Finished heliosFlutter.jl")