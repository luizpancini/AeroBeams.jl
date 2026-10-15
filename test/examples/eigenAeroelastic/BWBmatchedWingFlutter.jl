using AeroBeams

# Airspeed range [m/s]
URange = collect(60:1:180)

# Flight altitude
h = feet2meters(20e3)

# Aerodynamic solver
aeroSolver = Indicial()

# Flag to update airfoil parameters
updateAirfoilParameters = true

# Number of modes to track
nModes = 5

# Trim system solver 
NR = create_NewtonRaphson(trimSolverMethod=:LM,displayStatus=false)

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 60

# Initialize outputs
eigenproblem = Array{EigenProblem}(undef,length(URange))
eigensols = Array{EigenSolution}(undef,length(URange))
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))

# Initial guess solutions
x0Trim = zeros(0)
x0Eigen = zeros(0)

# Sweep airspeed range
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Model for trim problem
    BWBtrim = first(create_BWB(aeroSolver=aeroSolver,altitude=h,airspeed=U,δElevIsTrimVariable=true,thrustIsTrimVariable=true,updateAirfoilParameters=updateAirfoilParameters))
    # Solve the trim problem
    trimProblem = create_TrimProblem(model=BWBtrim,x0=x0Trim,systemSolver=NR)
    solve!(trimProblem)
    # Trim outputs
    trimSol = trimProblem.trimSolution
    aeroVars = trimProblem.loadstepHistory.aeroVariables
    trimAoA[i] = aeroVars[end][BWBtrim.beams[3].elements[1].assemblyData.globalID].unsteadyKin.αₑ
    trimThrust[i] = trimSol.trimLoads[1]
    trimδ[i] = trimSol.trimδ[1]
    println("Trim AoA = $(round(rad2deg(trimAoA[i]),digits=2)) deg, trim thrust = $(round(trimThrust[i],digits=2)) N, trim δ = $(round(rad2deg(trimδ[i]),digits=2)) deg")
    # Wing model for eigenproblem
    _,rightWingModel = create_BWB(aeroSolver=aeroSolver,altitude=h,airspeed=U,δElev=trimδ[i],wingRootAoA=trimAoA[i],thrust=trimThrust[i],updateAirfoilParameters=updateAirfoilParameters)
    # Solve eigenproblem
    eigenproblem[i] = create_EigenProblem(model=rightWingModel,eigenmethod=eigenmethod,nev=nev)
    solve!(eigenproblem[i])
    # Eigensolution
    eigensols[i] = eigenproblem[i].eigensolution
    # Set initial guess solutions for next iteration
    global x0Trim = trimProblem.systemArrays.x
    global x0Eigen = eigenproblem[i].systemArrays.x
end

# Apply mode tracking
freqs,damps,modeShapes = mode_tracking(eigensols)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,modeShapesOsc,_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=1e0, maxTracks=nModes, dropOutsideBounds=true)

# Flutter onset speeds, frequencies and modes
flutterSpeeds,flutterFreqs,flutterModes,flutterInds,_ = find_flutter(URange,freqsOsc,dampsOsc,1:length(URange))
for (mode,speed,freq) in zip(flutterModes,flutterSpeeds,flutterFreqs)
    println("Mode $mode: flutter speed = $(round(speed,sigdigits=4)) m/s, flutter frequency = $(round(rads2Hz(freq),digits=2)) Hz")
end

## Plots
# Set paths
relPath = "/test/outputs/figures/BWBmatchedWingFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Each mode shape animation at lowest flutter speed
iFlutter = round(Int,flutterInds[1])+1
for mode in 1:nModes
    modeAnim = plot_mode_shapes_animation(eigenproblem[iFlutter],modeShapes=modeShapesOsc[mode:mode,iFlutter],modeNumbers=[mode],nFramesPerCycle=21,fps=20,scale=2,view=(45,30),element2centralize=11,legendPos=:top,plotBCs=false,plotAxes=false,displayProgress=true,save=true,savePath=string(relPath,"/BWBmatchedWingFlutter_mode",mode,".gif"))
    display(modeAnim)
end

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
msw = 0
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=extrema(URange), tickfont=font(ts), guidefont=font(fs))
plot!(URange, rad2deg.(trimAoA), c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/BWBmatchedWingFlutter_AoA.pdf"))

# Trim thrust force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=extrema(URange), tickfont=font(ts), guidefont=font(fs))
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/BWBmatchedWingFlutter_thrust.pdf"))

# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=extrema(URange), tickfont=font(ts), guidefont=font(fs))
plot!(URange, rad2deg.(trimδ), c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/BWBmatchedWingFlutter_delta.pdf"))

# Root locus
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-30,5],ylims=[0,125], tickfont=font(ts), guidefont=font(fs))
for mode in 1:nModes
    scatter!(dampsOsc[mode,:], freqsOsc[mode,:], c=modeColors[mode], ms=ms, msw=msw, label=false)
end
display(plt_RL)
savefig(string(absPath,"/BWBmatchedWingFlutter_rootlocus.pdf"))

# V-g-f
plt51 = plot(ylabel="Frequency [rad/s]", xlims=extrema(URange), ylims=[0,100], tickfont=font(ts), guidefont=font(12))
for mode in 1:nModes
    plot!(URange, freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt52 = plot(xlabel="Airspeed [m/s]", ylabel="Damping ratio", xlims=extrema(URange), ylims=[-0.3,0.2], tickfont=font(ts), guidefont=font(12), yticks=-0.3:0.1:0.2)
for mode in 1:nModes
    plot!(URange, dampsOsc[mode,:]./freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt5 = plot(plt51,plt52, layout=(2,1))
display(plt5)
savefig(string(absPath,"/BWBmatchedWingFlutter_Vgf.pdf"))

println("Finished BWBmatchedWingFlutter.jl")