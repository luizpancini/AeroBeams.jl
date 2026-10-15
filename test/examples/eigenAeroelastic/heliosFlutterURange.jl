using AeroBeams

# Airspeed range [m/s]
URange = collect(5:0.5:20)

# Stiffness factor
λ = 1

# Payload [lb]
P = 100

# Flag to include beam pods
beamPods = true

# Flag to set payload on wing
payloadOnWing = false

# Aerodynamic solver
aeroSolver = Indicial()

# Attachment spring stiffness
μ = 1e-1

# System solver for trim problem
NR = create_NewtonRaphson(trimSolverMethod=:LQ,displayStatus=false)

# Number of modes to track
nModes = 10

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :arnoldi
nev = 150

# Initialize outputs
eigenProblem = Array{EigenProblem}(undef,length(URange))
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))
eigensols = Array{EigenSolution}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Model for trim problem
    heliosTrim,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,payloadOnWing=payloadOnWing,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true,withAttachmentSprings=true,μ=μ)
    # Create and solve trim problem
    trimProblem = create_TrimProblem(model=heliosTrim,systemSolver=NR)
    solve!(trimProblem)
    # Extract trim variables
    trimSol = trimProblem.trimSolution
    aeroVars = trimProblem.loadstepHistory.aeroVariables
    trimAoA[i] = aeroVars[end][midSpanElem].unsteadyKin.αₑ
    trimThrust[i] = trimSol.trimLoads[1]
    trimδ[i] = trimSol.trimδ[1]
    println("Trim AoA = $(round(rad2deg(trimAoA[i]), digits=2)) deg, trim thrust = $(round(trimThrust[i], digits=2)) N, trim δ = $(round(rad2deg(trimδ[i]), digits=2)) deg")
    # Model for eigen problem
    heliosEigen,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,payloadOnWing=payloadOnWing,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δ=trimδ[i],thrust=trimThrust[i],withAttachmentSprings=true,μ=μ)
    # Set initial guess solution
    x0 = i==1 ? zeros(0) : eigenProblem[i-1].systemArrays.x
    # Create and solve eigen problem
    eigenProblem[i] = create_EigenProblem(model=heliosEigen,eigenmethod=eigenmethod,nev=nev,refTrimProblem=trimProblem)
    solve_eigen!(eigenProblem[i])
    # Eigensolution
    eigensols[i] = eigenProblem[i].eigensolution
end

# Apply mode tracking
freqs,damps,modeShapes,_ = mode_tracking(eigensols)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,modeShapesOsc,_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=1e-1, maxTracks=nModes)


# Bending and torsional mode shapes at selected airspeed
ind = findfirst(U -> U >= 15, URange)
nElemWing = sum([eigenProblem[1].model.beams[i].discretization.nElements for i in 1:4])
u2_modeShapes = [reduce(vcat,
                 [modeShapesOsc[m,ind].state[e].nodal.u1[2],
                  modeShapesOsc[m,ind].state[e].nodal.u2[2]] 
                  for e in 1:nElemWing) for m in 1:nModes]
u3_modeShapes = [reduce(vcat,
                 [modeShapesOsc[m,ind].state[e].nodal.u1[3],
                  modeShapesOsc[m,ind].state[e].nodal.u2[3]] 
                  for e in 1:nElemWing) for m in 1:nModes]
p1_modeShapes = [reduce(vcat,
                 [modeShapesOsc[m,ind].state[e].nodal.p1[1],
                  modeShapesOsc[m,ind].state[e].nodal.p2[1]] 
                  for e in 1:nElemWing) for m in 1:nModes]

## Plots
# Set paths
relPath = "/test/outputs/figures/heliosFlutterURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at the last airspeed
modesPlot = plot_mode_shapes(eigenProblem[end],modeShapes=modeShapesOsc[:,end],modeNumbers=collect(1:nModes),scale=5,view=(30,30),legendPos=(0.3,0.98),save=true,savePath=string(relPath,"/heliosFlutterURange_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [rad/s]", xlims=extrema(URange))
for mode in 1:nModes
    plot!(URange, freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URange), ylims=[-1.0,0.25], legend=:bottomleft)
for mode in 1:nModes
    plot!(URange, dampsOsc[mode,:]./freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/heliosFlutterURange_Vgf.pdf"))

# Root locus
plt2 = plot(xlabel="Damping ratio", ylabel="Frequency [rad/s]", xlims=[-1.0,0.25])
for mode in 1:nModes
    plot!(dampsOsc[mode,:]./freqsOsc[mode,:], freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
display(plt2)
savefig(string(absPath,"/heliosFlutterURange_rootlocus.pdf"))

# Arclength of wing
x1 = Float64[]
for i=1:4
    beam = eigenProblem[1].model.beams[i]
    nElem = beam.discretization.nElements
    for e in 1:nElem
        append!(x1, beam.elements[e].geometry.x1_n1, beam.elements[e].geometry.x1_n2)
    end
    if i > 1
        x1[end-2*nElem+1:end] .+= x1[end-2*nElem]
    end
end
L = maximum(x1)

# Mode shapes at specified velocity
plt4 = plot(xlabel="\$x1/L\$", ylabel="\$u_2\$", title="Chordwise bending mode shape at U=$(URange[ind]) m/s")
plt5 = plot(xlabel="\$x1/L\$", ylabel="\$u_3\$", title="Flapwise bending mode shape at U=$(URange[ind]) m/s")
plt6 = plot(xlabel="\$x1/L\$", ylabel="\$p_1\$", title="Torsional mode shape at U=$(URange[ind]) m/s")
for m in 1:nModes
    plot!(plt4, x1/L, u2_modeShapes[m], c=modeColors[m], lw=lw, label=false)
    plot!(plt5, x1/L, u3_modeShapes[m], c=modeColors[m], lw=lw, label=false)
    plot!(plt6, x1/L, p1_modeShapes[m], c=modeColors[m], lw=lw, label=false)
end
display(plt4)
display(plt5)
display(plt6)
savefig(string(absPath,"/SMWFlutter_modalDisp.pdf"))

println("Finished heliosFlutterURange.jl")