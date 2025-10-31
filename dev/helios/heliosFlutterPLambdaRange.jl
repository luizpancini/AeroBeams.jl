using AeroBeams

# Option for reduced chord
reducedChord = false

# Option to include beam pods
beamPods = true

# Option to set payload on wing
payloadOnWing = false

# Circulatory indicial function
circulatoryIndicialFunction = "Wagner"

# Aerodynamic solver
aeroSolver = Indicial(circulatoryIndicialFunction=circulatoryIndicialFunction)

# Wing airfoil
wingAirfoil = deepcopy(NACA23012A)

# Stiffness range
λRange = [1, 1e4]

# Payload range
PRange = collect(0:20:500)

# Airspeed
U = 40*0.3048

# Discretization
nElemStraightSemispan = 10
nElemDihedralSemispan = 5
nElemPod = 2

# Value for "rigid" structural properties (not rigid aircraft properties)
∞ = 1e9

# Number of modes
nModes = 10

# System solver for trim problem
relaxFactor = 0.5
maxIter = 100
relTol = 1e-12
absTol = 1e-12
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,relativeTolerance=relTol,absoluteTolerance=absTol)

# Attachment springs
μ = 1e-3
ku = μ*[1; 1; 1]
kp = ku
spring = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)

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

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Sweep payload
    for (j,P) in enumerate(PRange)
        # Display progress
        println("Solving for λ = $λ, payload = $P lb")
        # Model for trim problem
        heliosTrim,midSpanElem,_,_,rightWingStraight,_ = create_Helios(aeroSolver=aeroSolver,wingAirfoil=wingAirfoil,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,reducedChord=reducedChord,payloadOnWing=payloadOnWing,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,∞=∞)
        # Add springs at wing root
        add_springs_to_beam!(beam=rightWingStraight,springs=[spring])
        # Update model
        update_model!(heliosTrim)
        # Set initial guess solution as previous known solution
        x0Trim = (j==1) ? zeros(0) : trimProblem[i,j-1].x
        # Create and solve trim problem
        trimProblem[i,j] = create_TrimProblem(model=heliosTrim,systemSolver=NRtrim,x0=x0Trim)
        solve!(trimProblem[i,j])
        # Extract trim variables
        trimAoA = trimProblem[i,j].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.αₑ*180/π
        trimThrust = trimProblem[i,j].x[end-1]*trimProblem[i,j].model.forceScaling
        trimδ = trimProblem[i,j].x[end]
        println("AoA = $(trimAoA), T = $(trimThrust), δ = $(trimδ*180/π)")
        # Model for eigen problem
        heliosEigen,_,_,_,rightWingStraight,_ = create_Helios(aeroSolver=aeroSolver,wingAirfoil=wingAirfoil,beamPods=beamPods,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δ=trimδ,thrust=trimThrust,reducedChord=reducedChord,payloadOnWing=payloadOnWing,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan,nElemPod=nElemPod,∞=∞)
        # Add springs at wing root
        add_springs_to_beam!(beam=rightWingStraight,springs=[spring])
        # Update model
        update_model!(heliosEigen)
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=heliosEigen,nModes=nModes,frequencyFilterLimits=[2e-2,Inf],refTrimProblem=trimProblem[i,j])
        solve_eigen!(eigenProblem[i,j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = eigenProblem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = eigenProblem[i,j].dampingsOscillatory
        untrackedEigenvectors[i,j] = eigenProblem[i,j].eigenvectorsOscillatoryCplx
    end
    # Apply mode tracking
    freqs[i,:],damps[i,:],_,matchedModes = mode_tracking_hungarian(PRange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
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
DPI = 300
fps = 5
labels = ["Flexible" "Rigid"]
mshape = [:circle, :star]
grad_flex = :cool
grad_rigid = :copper
colors_flex = cgrad(grad_flex, length(PRange), categorical=true)
colors_rigid = cgrad(grad_rigid, length(PRange), categorical=true)
if typeof(aeroSolver) == Indicial
    aeroSolverLabel = "_linear"
elseif typeof(aeroSolver) == BLi
    aeroSolverLabel = "_ds"
end

# Root locus 
p1 = scatter([NaN], [NaN], zcolor=[NaN], clims=(0,500), label=false, c=grad_flex, colorbar_titlefontsize=12, colorbar_title="\nPayload [lb] - Flexible", right_margin=5.5Plots.mm, xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,1], ylims=[0,10], tickfont=font(ts), guidefont=font(fs))
p2 = scatter([NaN], [NaN], zcolor=[NaN], clims=(0,500), framestyle=:none, axis=false, label=false, grid=false, c=grad_rigid, colorbar_titlefontsize=12, colorbar_title="\n Payload [lb] - Rigid", right_margin=5.5Plots.mm, tickfont=font(ts))
l = @layout [grid(1,1) a{0.01w}]
plt_RL = plot(p1,p2,layout=l)
plot!([0; 0], [0, 10], c=:black, ls=:dash, lw=2, label=false)
for mode in 1:nModes
    for (j,P) in enumerate(PRange)
        scatter!([modeDampings[1,mode][j]], [modeFrequencies[1,mode][j]], c=colors_flex[j], shape=mshape[1], ms=2, msw=msw, label=false)
        scatter!([modeDampings[2,mode][j]], [modeFrequencies[2,mode][j]], c=colors_rigid[j], shape=mshape[2], ms=2, msw=msw, label=false)
    end
end
display(plt_RL)
savefig(plt_RL,string(absPath,"/heliosFlutterPLambdaRange_rootlocus",aeroSolverLabel,".pdf"))

# Animated root locus plot
plot_RL_anim = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,1], ylims=[0,10], tickfont=font(ts), guidefont=font(fs), dpi=DPI)
plot!([0; 0], [0, 10], c=:black, ls=:dash, lw=2, label=false)
scatter!([NaN], [NaN], c=colors_flex[1], shape=mshape[1], ms=ms, msw=msw, label=labels[1])
scatter!([NaN], [NaN], c=colors_rigid[1], shape=mshape[2], ms=ms, msw=msw, label=labels[2])
anim = @animate for (j,P) in enumerate(PRange)
    title!("\$P = $P\$ lb")
    for mode in 1:nModes
        plot!([modeDampings[1,mode][j]], [modeFrequencies[1,mode][j]], c=colors_flex[j], shape=mshape[1], ms=ms, msw=msw, label=false)
        plot!([modeDampings[2,mode][j]], [modeFrequencies[2,mode][j]], c=colors_rigid[j], shape=mshape[2], ms=ms, msw=msw, label=false)
    end
end
gif_handle = gif(anim, string(absPath,"/heliosFlutterPLambdaRange_rootlocus.gif"), fps=fps)
display(gif_handle)

# Root locus (phugoid zoom)
plt_RLphugoid = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.05,0.2], ylims=[0,0.8], tickfont=font(ts), guidefont=font(fs))
plot!([0; 0], [0, 10], c=:black, ls=:dash, lw=2, label=false)
for mode in 1:nModes
    for (j,P) in enumerate(PRange)
        scatter!([modeDampings[1,mode][j]],[ modeFrequencies[1,mode][j]], c=colors_flex[j], shape=mshape[1], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[2,mode][j]], [modeFrequencies[2,mode][j]], c=colors_rigid[j], shape=mshape[2], ms=ms, msw=msw, label=false)
    end
end
display(plt_RLphugoid)
savefig(plt_RLphugoid,string(absPath,"/heliosFlutterPLambdaRange_phugoid",aeroSolverLabel,".pdf"))

# Animated root locus (zoom on phugoid)
plot_RL_anim = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.05,0.2], ylims=[0,0.8], tickfont=font(ts), guidefont=font(fs), dpi=DPI)
plot!([0; 0], [0, 10], c=:black, ls=:dash, lw=2, label=false)
anim = @animate for (j,P) in enumerate(PRange)
    title!("\$P = $P\$ lb")
    for mode in 1:nModes
        plot!([modeDampings[1,mode][j]], [modeFrequencies[1,mode][j]], c=colors_flex[j], shape=mshape[1], ms=ms, msw=msw, label=false)
        plot!([modeDampings[2,mode][j]], [modeFrequencies[2,mode][j]], c=colors_rigid[j], shape=mshape[2], ms=ms, msw=msw, label=false)
    end
end
gif_handle = gif(anim, string(absPath,"/heliosFlutterPLambdaRange_rootlocus_phugoid.gif"), fps=fps)
display(gif_handle)

# Root locus (zoom on low frequency modes)
plt_RLlow = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.75,0.2], ylims=[0,0.8], tickfont=font(ts), guidefont=font(fs))
plot!([0; 0], [0, 10], c=:black, ls=:dash, lw=2, label=false)
for mode in 1:nModes
    for (j,P) in enumerate(PRange)
        scatter!([modeDampings[1,mode][j]],[ modeFrequencies[1,mode][j]], c=colors_flex[j], shape=mshape[1], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[2,mode][j]], [modeFrequencies[2,mode][j]], c=colors_rigid[j], shape=mshape[2], ms=ms, msw=msw, label=false)
    end
end
display(plt_RLlow)
savefig(plt_RLlow,string(absPath,"/heliosFlutterPLambdaRange_low",aeroSolverLabel,".pdf"))

# Root locus (zoom on short period)
plt_RLsp = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,-2.5], ylims=[2,4], tickfont=font(ts), guidefont=font(fs))
plot!([0; 0], [0, 10], c=:black, ls=:dash, lw=2, label=false)
for mode in 1:nModes
    for (j,P) in enumerate(PRange)
        scatter!([modeDampings[1,mode][j]],[ modeFrequencies[1,mode][j]], c=colors_flex[j], shape=mshape[1], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[2,mode][j]], [modeFrequencies[2,mode][j]], c=colors_rigid[j], shape=mshape[2], ms=ms, msw=msw, label=false)
    end
end
display(plt_RLsp)
savefig(plt_RLsp,string(absPath,"/heliosFlutterPLambdaRange_sp",aeroSolverLabel,".pdf"))

# Mode shapes of flexible aircraft at highest payload
modesPlot_flex = plot_mode_shapes(eigenProblem[1,end],element2centralize=nElemStraightSemispan+nElemDihedralSemispan,scale=10,view=(30,30),legendPos=:topright,nModes=6,modalColorScheme=:rainbow,save=true,savePath=string(relPath,"/heliosFlutterPLambdaRange_modeShapesFlex",aeroSolverLabel,".pdf"))
display(modesPlot_flex)

# Static mode shapes of rigid aircraft at highest payload
modesPlot_rigid = plot_mode_shapes(eigenProblem[2,end],element2centralize=nElemStraightSemispan+nElemDihedralSemispan,scale=10,view=(45,15),legendPos=:topright,nModes=4,modeLabels=["Dutch roll","Phugoid","Roll","Short-period"],modalColorScheme=:rainbow,save=true,savePath=string(relPath,"/heliosFlutterPLambdaRange_modeShapesRigid",aeroSolverLabel,".pdf"))
display(modesPlot_rigid)

# Mode shapes animation of flexible aircraft at selected payload
payloadAnim = 200
ind = findfirst(x -> x >= payloadAnim, PRange)
modesAnim_flex = plot_mode_shapes_animation(eigenProblem[1,ind],element2centralize=nElemStraightSemispan+nElemDihedralSemispan,matchModeFrequency=false,scale=20,view=(45,15),nFramesPerCycle=51,displayProgress=true,plotBCs=false,plotAxes=false,showLegend=false,legendPos=:top,save=true,savePath=string(relPath,"/heliosFlutterPLambdaRange_modeShapesFlex",aeroSolverLabel,".gif"))
display(modesAnim_flex)
modes2Anim = [3]
for (m,mode) in enumerate(modes2Anim)
    modeAnim_flex = plot_mode_shapes_animation(eigenProblem[1,ind],element2centralize=nElemStraightSemispan+nElemDihedralSemispan,modes2plot=[mode],numberOfColors=nModes,scale=20,view=(45,15),nFramesPerCycle=21,displayProgress=true,plotBCs=false,plotAxes=false,showLegend=false,legendPos=:top,save=true,savePath=string(relPath,"/heliosFlutterPLambdaRange_modeShape",mode,"Flex",aeroSolverLabel,".gif"))
    display(modeAnim_flex)
end

# # Mode shapes animation of rigid aircraft at highest payload
# modesAnim_rigid = plot_mode_shapes_animation(eigenProblem[2,end],element2centralize=nElemStraightSemispan+nElemDihedralSemispan,scale=20,view=(45,15),nModes=4,numberOfColors=4,nFramesPerCycle=101,displayProgress=true,plotBCs=false,plotAxes=false,legendFontSize=10,legendPos=:top,modeLabels=["Dutch roll","Phugoid","Roll","Short-period"],save=true,savePath=string(relPath,"/heliosFlutterPLambdaRange_modeShapesRigid",aeroSolverLabel,".gif"))

println("Finished heliosFlutterPLambdaRange.jl")