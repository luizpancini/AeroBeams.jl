using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Mode tracking option
modeTracking = true

# Stiffness factor
λ = 1e0

# Aerodynamic solver
aeroSolver = Indicial()

# Tip loss option
hasTipCorrection = true

# Option to update airfoil parameters (apply Prandtl-Glauert lift-slope correction)
updateAirfoilParameters = false

# Set NR system solver 
relaxFactor = 0.5
displayStatus = false
maxiter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxiter,displayStatus=displayStatus)

# Set number of vibration modes
nModes = 8

# Set airspeed range and initialize outputs
URange = vcat(30:2:160)
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))

# Attachment springs
μ = 1e-2
ku = μ*[1; 1; 1]
kp = ku
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)
spring2 = create_Spring(elementsIDs=[3],nodesSides=[2],ku=ku,kp=kp)

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Model for trim problem
    BWBtrim = create_BWB(aeroSolver=aeroSolver,stiffnessFactor=λ,hasTipCorrection=hasTipCorrection,updateAirfoilParameters=updateAirfoilParameters,airspeed=U,δElevIsTrimVariable=true,thrustIsTrimVariable=true)
    # Add springs
    add_springs_to_beam!(beam=BWBtrim.beams[2],springs=[spring1])
    add_springs_to_beam!(beam=BWBtrim.beams[3],springs=[spring2])
    # Update model
    BWBtrim.skipValidationMotionBasisA = true
    update_model!(BWBtrim)
    # Set initial guess solution as previous known solution
    x0Trim = i == 1 ? zeros(0) : trimProblem.x
    # Create and trim problem
    global trimProblem = create_TrimProblem(model=BWBtrim,systemSolver=NR,x0=x0Trim)
    solve!(trimProblem)
    # Extract trim variables
    trimAoA[i] = trimProblem.aeroVariablesOverσ[end][BWBtrim.beams[3].elementRange[1]].flowAnglesAndRates.αₑ
    trimThrust[i] = trimProblem.x[end-1]*BWBtrim.forceScaling 
    trimδ[i] = trimProblem.x[end]
    println("Trim AoA = $(trimAoA[i]*180/π), trim thrust = $(trimThrust[i]), trim δ = $(trimδ[i]*180/π)")
    # Model for eigen problem
    BWBeigen = create_BWB(aeroSolver=aeroSolver,stiffnessFactor=λ,hasTipCorrection=hasTipCorrection,updateAirfoilParameters=updateAirfoilParameters,airspeed=U,δElev=trimδ[i],thrust=trimThrust[i])
    # Add springs
    add_springs_to_beam!(beam=BWBeigen.beams[2],springs=[spring1])
    add_springs_to_beam!(beam=BWBeigen.beams[3],springs=[spring2])
    # Update model
    BWBeigen.skipValidationMotionBasisA = true
    update_model!(BWBeigen)
    # Create and solve eigen problem
    global eigenProblem = create_EigenProblem(model=BWBeigen,nModes=nModes,frequencyFilterLimits=[1e0,Inf64],jacobian=trimProblem.jacobian[1:end,1:end-trimProblem.model.nTrimVariables],inertia=trimProblem.inertia)
    solve_eigen!(eigenProblem)
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = eigenProblem.frequenciesOscillatory
    untrackedDamps[i] = round_off!(eigenProblem.dampingsOscillatory,1e-12)
    untrackedEigenvectors[i] = eigenProblem.eigenvectorsOscillatoryCplx
end

# Frequencies and dampings after mode tracking
if modeTracking
    freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
else
    freqs,damps = untrackedFreqs,untrackedDamps
end

# Separate frequencies and damping ratios by mode
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeFrequencies =  Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
end

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(URange)))
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
relPath = "/test/outputs/figures/BWBstability"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Mode shapes at highest airspeed
modesPlot = plot_mode_shapes(eigenProblem,scale=1,view=(30,30),legendPos=:outerright,save=true,savePath=string(relPath,"/BWBstability_modeShapes.pdf"))
display(modesPlot)
# Trim root angle of attack vs airspeed
gr()
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[URange[1],URange[end]])
plot!(URange, trimAoA*180/π, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/BWBstability_AoA.pdf"))
# Trim propeller force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[URange[1],URange[end]])
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/BWBstability_thrust.pdf"))
# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[URange[1],URange[end]])
plot!(URange, trimδ*180/π, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/BWBstability_delta.pdf"))
# Root locus
plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-20,5],ylims=[0,120])
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
display(plt4)
savefig(string(absPath,"/BWBstability_rootlocus.pdf"))
# V-g-f
plt51 = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,120])
for mode in 1:nModes
    scatter!(URange, modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=0,  label=false)
end
plt52 = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", xlims=[URange[1],URange[end]], ylims=[-10,5])
for mode in 1:nModes
    scatter!(URange, modeDampings[mode], c=modeColors[mode], ms=ms, msw=0, label=false)
end
plt5 = plot(plt51,plt52, layout=(2,1))
display(plt5)
savefig(string(absPath,"/BWBstability_Vgf.pdf"))

println("Finished BWBstability.jl")