using AeroBeams

# Option for reduced chord
reducedChord = true

# Option to include beam pods
beamPods = true

# Option to set payload on wing
payloadOnWing = false

# Aerodynamic solver
aeroSolver = Indicial()

# Wing airfoil
wingAirfoil = deepcopy(NACA23012A)

# Stiffness range
λRange = [1, 1e4]

# Payload range
PRange = collect(0:10:500)

# Airspeed
U = 40*0.3048

# Option for mode tracking
modeTracking = true

# Number of modes
nModes = 10

# System solver for trim problem
relaxFactor = 0.5
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=100,displayStatus=false)

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(λRange),length(PRange))
eigenProblem = Array{EigenProblem}(undef,length(λRange),length(PRange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(λRange),length(PRange))
untrackedDamps = Array{Vector{Float64}}(undef,length(λRange),length(PRange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(λRange),length(PRange))
freqs = Array{Vector{Float64}}(undef,length(λRange),length(PRange))
damps = Array{Vector{Float64}}(undef,length(λRange),length(PRange))
modeFrequencies =  Array{Vector{Float64}}(undef,length(λRange),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(λRange),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(λRange),nModes)

# Attachment springs
μ = 1e-2
ku = μ*[1; 1; 1]
kp = ku
spring = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Sweep payload
    for (j,P) in enumerate(PRange)
        # Display progress
        println("Solving for λ = $λ, payload = $P lb")
        # Model for trim problem
        heliosTrim,midSpanElem,_,_,rightWingStraight,_ = create_Helios(aeroSolver=aeroSolver,wingAirfoil=wingAirfoil,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,reducedChord=reducedChord,payloadOnWing=payloadOnWing)
        # Add springs at wing root
        add_springs_to_beam!(beam=rightWingStraight,springs=[spring])
        # Update model
        heliosTrim.skipValidationMotionBasisA = true
        update_model!(heliosTrim)
        # Set initial guess solution as previous known solution
        x0Trim = (j==1) ? zeros(0) : trimProblem[i,j-1].x
        # Create and solve trim problem
        trimProblem[i,j] = create_TrimProblem(model=heliosTrim,systemSolver=NR,x0=x0Trim)
        solve!(trimProblem[i,j])
        # Extract trim variables
        trimAoA = trimProblem[i,j].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.αₑ*180/π
        trimThrust = trimProblem[i,j].x[end-1]*trimProblem[i,j].model.forceScaling
        trimδ = trimProblem[i,j].x[end]
        println("AoA = $(trimAoA), T = $(trimThrust), δ = $(trimδ*180/π)")
        # Model for eigen problem
        heliosEigen,_,_,_,rightWingStraight,_ = create_Helios(aeroSolver=aeroSolver,wingAirfoil=wingAirfoil,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δ=trimδ,thrust=trimThrust,reducedChord=reducedChord,payloadOnWing=payloadOnWing)
        # Add springs at wing root
        add_springs_to_beam!(beam=rightWingStraight,springs=[spring])
        # Update model
        heliosEigen.skipValidationMotionBasisA = true
        update_model!(heliosEigen)
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=heliosEigen,nModes=nModes,frequencyFilterLimits=[2e-2,Inf64],jacobian=trimProblem[i,j].jacobian[1:end,1:end-2],inertia=trimProblem[i,j].inertia,refTrimProblem=trimProblem[i,j])
        solve_eigen!(eigenProblem[i,j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = eigenProblem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(eigenProblem[i,j].dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = eigenProblem[i,j].eigenvectorsOscillatoryCplx
    end
    # Apply mode tracking, if applicable
    if modeTracking
        freqs[i,:],damps[i,:],_,matchedModes = mode_tracking(PRange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    else
        freqs[i,:],damps[i,:] = untrackedFreqs[i,:],untrackedDamps[i,:]
    end
    # Separate frequencies and dampings by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(PRange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(PRange)]
        modeDampingRatios[i,mode] = modeDampings[i,mode]./modeFrequencies[i,mode]
    end
end

# Set paths
relPath = "/dev/helios/figures/heliosFlutterPLambdaRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
using Plots, ColorSchemes
gr()
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
msw = 0
labels = ["Flexible" "Rigid"]
mshape = [:circle, :star]
grad_flex = :cool
grad_rigid = :copper
colors_flex = cgrad(grad_flex, length(PRange), categorical=true)
colors_rigid = cgrad(grad_rigid, length(PRange), categorical=true)

# Root locus 
p1 = scatter([NaN], [NaN], zcolor=[NaN], clims=(0,500), label=false, c=grad_flex, colorbar_titlefontsize=12, colorbar_title="\nPayload [lb] - Flexible", right_margin=5.5Plots.mm, xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,1], ylims=[0,10], tickfont=font(ts), guidefont=font(fs))
p2 = scatter([NaN], [NaN], zcolor=[NaN], clims=(0,500), framestyle=:none, axis=false, label=false, grid=false, c=grad_rigid, colorbar_titlefontsize=12, colorbar_title="\n Payload [lb] - Rigid", right_margin=5.5Plots.mm, tickfont=font(ts))
l = @layout [grid(1,1) a{0.01w}]
plt_RL = plot(p1,p2,layout=l)
display(plt_RL)
for mode in 1:nModes
    for (j,P) in enumerate(PRange)
        scatter!([modeDampings[1,mode][j]], [modeFrequencies[1,mode][j]], c=colors_flex[j], shape=mshape[1], ms=2, msw=msw, label=false)
        scatter!([modeDampings[2,mode][j]], [modeFrequencies[2,mode][j]], c=colors_rigid[j], shape=mshape[2], ms=2, msw=msw, label=false)
    end
end
display(plt_RL)
savefig(string(absPath,"/heliosFlutterPLambdaRange_rootlocus.pdf"))

# Root locus (phugoid zoom)
plt_RLphugoid = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.05,0.2], ylims=[0,0.8], tickfont=font(ts), guidefont=font(fs))
for mode in 1:nModes
    for (j,P) in enumerate(PRange)
        scatter!([modeDampings[1,mode][j]],[ modeFrequencies[1,mode][j]], c=colors_flex[j], shape=mshape[1], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[2,mode][j]], [modeFrequencies[2,mode][j]], c=colors_rigid[j], shape=mshape[2], ms=ms, msw=msw, label=false)
    end
end
display(plt_RLphugoid)
savefig(string(absPath,"/heliosFlutterPLambdaRange_phugoid.pdf"))

# Root locus (zoom on low frequency modes)
plt_RLlow = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.75,0.2], ylims=[0,0.8], tickfont=font(ts), guidefont=font(fs))
for mode in 1:nModes
    for (j,P) in enumerate(PRange)
        scatter!([modeDampings[1,mode][j]],[ modeFrequencies[1,mode][j]], c=colors_flex[j], shape=mshape[1], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[2,mode][j]], [modeFrequencies[2,mode][j]], c=colors_rigid[j], shape=mshape[2], ms=ms, msw=msw, label=false)
    end
end
display(plt_RLlow)
savefig(string(absPath,"/heliosFlutterPLambdaRange_low.pdf"))

# Root locus (zoom on short period)
plt_RLsp = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,-2.5], ylims=[2,4], tickfont=font(ts), guidefont=font(fs))
for mode in 1:nModes
    for (j,P) in enumerate(PRange)
        scatter!([modeDampings[1,mode][j]],[ modeFrequencies[1,mode][j]], c=colors_flex[j], shape=mshape[1], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[2,mode][j]], [modeFrequencies[2,mode][j]], c=colors_rigid[j], shape=mshape[2], ms=ms, msw=msw, label=false)
    end
end
display(plt_RLsp)
savefig(string(absPath,"/heliosFlutterPLambdaRange_sp.pdf"))

# Mode shapes of flexible aircraft at highest payload
modesPlot_flex = plot_mode_shapes(eigenProblem[1,end],interactive=false,ΔuDef=-eigenProblem[1,end].elementalStatesOverσ[end][15].u,scale=10,view=(30,30),legendPos=:topright,nModes=6,modalColorScheme=:rainbow,save=true,savePath=string(relPath,"/heliosFlutterPLambdaRange_modeShapesFlex.pdf"))
display(modesPlot_flex)

# Mode shapes of rigid aircraft at highest payload
modesPlot_rigid = plot_mode_shapes(eigenProblem[2,end],interactive=false,ΔuDef=-eigenProblem[2,end].elementalStatesOverσ[end][15].u,scale=10,view=(45,15),legendPos=:topright,nModes=5,modeLabels=["Dutch roll","Phugoid","Roll","Short-period","1st OOP"],modalColorScheme=:rainbow,save=true,savePath=string(relPath,"/heliosFlutterPLambdaRange_modeShapesRigid.pdf"))
display(modesPlot_rigid)

println("Finished heliosFlutterPLambdaRange.jl")