using AeroBeams

# Airspeed range [m/s]
URange = vcat(1e-3,0.5:0.5:50)

# Pitch angle [rad]
θ = deg2rad(0)

# Altitude [m]
h = 20e3

# Discretization
nElem = 32

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil
airfoil = deepcopy(FlatPlate)

# Model
SMWFlutter,L = create_SMW(aeroSolver=aeroSolver,airfoil=airfoil,pitchAngle=θ,altitude=h,nElem=nElem)

# Undeformed nodal and midpoint positions
wing = SMWFlutter.beams[1]
x1_0 = vcat([vcat(wing.elements[e].assemblyData.r_n1[1],wing.elements[e].assemblyData.r_n2[1]) for e in 1:nElem]...)
x3_0 = -vcat([vcat(wing.elements[e].assemblyData.r_n1[3],wing.elements[e].assemblyData.r_n2[3]) for e in 1:nElem]...)
x1_e = [wing.elements[e].geometry.x1 for e in 1:nElem]

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 80

# Initialize outputs
problem = Array{EigenProblem}(undef,length(URange))
eigensols = Array{EigenSolution}(undef,length(URange))
α_of_x1 = Array{Vector{Float64}}(undef,length(URange))
tip_u3 = Array{Float64}(undef,length(URange))
x1_def = Array{Vector{Float64}}(undef,length(URange))
x3_def = Array{Vector{Float64}}(undef,length(URange))

# Initialize guess solution
x0 = zeros(0)

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update airspeed
    update_v_A!(SMWFlutter,v_A=[0,U,0])
    # Create and solve problem
    problem[i] = create_EigenProblem(model=SMWFlutter,eigenmethod=eigenmethod,nev=nev,x0=x0)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
    # Final states
    state = problem[i].loadstepHistory.states[end]
    aeroVars = problem[i].loadstepHistory.aeroVariables[end]
    # Angle of attack over span
    α_of_x1[i] = [rad2deg(aeroVars[e].unsteadyKin.αₑ) for e in 1:nElem]
    # Tip OOP displacement
    tip_u3[i] = state[nElem].nodal.u2_b[3]
    # Displacements over span
    u1_of_x1 = vcat([vcat(state[e].nodal.u1_b[1],state[e].nodal.u2_b[1]) for e in 1:nElem]...)
    u2_of_x1 = vcat([vcat(state[e].nodal.u1_b[2],state[e].nodal.u2_b[2]) for e in 1:nElem]...)
    u3_of_x1 = vcat([vcat(state[e].nodal.u1_b[3],state[e].nodal.u2_b[3]) for e in 1:nElem]...)
    # Deformed nodal positions
    x1_def[i] = x1_0 .+ u1_of_x1
    x3_def[i] = x3_0 .+ u3_of_x1
    # Set initial guess solution for next iteration
    global x0 = problem[i].systemArrays.x
end

# Apply mode tracking
freqs,damps,modeShapes,_ = mode_tracking(eigensols)

# Get frequencies and dampings of oscillatory and non-oscillatory modes
freqsOsc,dampsOsc,modeShapesOsc,_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=1e-2, maxTracks=nModes)
freqsNonOsc,dampsNonOsc,_ = filter_tracked_modes(freqs,damps; fmin=0, fmax=0)

# Divergence speed
divergenceSpeeds,_ = find_divergence(URange,freqsNonOsc,dampsNonOsc)
divergenceSpeed = minimum(divergenceSpeeds)
println("Divergence speed = $(round(divergenceSpeed, digits=2)) m/s")

# Flutter onset and offset speeds, frequencies, and modes
flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetModes,flutterOnsetInds,flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetModes,flutterOffsetInds = find_flutter(URange,freqsOsc,dampsOsc,1:length(URange))

# Show flutter results
for (m,speed,freq) in zip(flutterOnsetModes,flutterOnsetSpeeds,flutterOnsetFreqs)
    println("Mode $m: onset speed = $(round(speed, digits=2)) m/s, frequency = $(round(freq, digits=2)) rad/s")
end
for (m,speed,freq) in zip(flutterOffsetModes,flutterOffsetSpeeds,flutterOffsetFreqs)
    println("Mode $m: offset speed = $(round(speed, digits=2)) m/s, frequency = $(round(freq, digits=2)) rad/s")
end

# Bending and torsional mode shapes at first flutter condition
iFlutter = round(Int, flutterOnsetInds[1])+1
u2_modeShapes = [reduce(vcat,
                 [modeShapesOsc[m,iFlutter].state[e].nodal.u1[2],
                  modeShapesOsc[m,iFlutter].state[e].nodal.u2[2]] 
                  for e in 1:nElem) for m in 1:nModes]
u3_modeShapes = [reduce(vcat,
                 [modeShapesOsc[m,iFlutter].state[e].nodal.u1[3],
                  modeShapesOsc[m,iFlutter].state[e].nodal.u2[3]] 
                  for e in 1:nElem) for m in 1:nModes]
p1_modeShapes = [reduce(vcat,
                 [modeShapesOsc[m,iFlutter].state[e].nodal.p1[1],
                  modeShapesOsc[m,iFlutter].state[e].nodal.p2[1]] 
                  for e in 1:nElem) for m in 1:nModes]

## Plots
# Set paths
relPath = "/test/outputs/figures/SMWFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at first flutter condition
modesPlot = plot_mode_shapes(problem[iFlutter],scale=5,modeShapes=modeShapesOsc[:,iFlutter],modeNumbers=1:nModes,view=(30,30),legendPos=(0.3,0.95),save=true,savePath=string(relPath,"/SMWFlutter_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
gr()

# Normalized deformed wingspan
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$", xlims=[0,1])
for (i,U) in enumerate(URange)
    plot!(x1_def[i]/L, x3_def[i]/L, lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt1)
savefig(string(absPath,"/SMWFlutter_disp.pdf"))

# Angle of attack over wingspan
plt11 = plot(xlabel="\$x_1/L\$", ylabel="\$\\alpha\$ [deg]", xlims=[0,1])
for (i,U) in enumerate(URange)
    plot!(x1_e/L, α_of_x1[i], lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt11)
savefig(string(absPath,"/SMWFlutter_AoA.pdf"))

# Root locus
plt2 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,1], ylims=[0,50])
for mode in 1:nModes
    plot!(dampsOsc[mode,:], freqsOsc[mode,:], c=modeColors[mode], lw=lw, label="Mode $mode")
end
display(plt2)
savefig(string(absPath,"/SMWFlutter_rootlocus.pdf"))

# V-g-f
plt31 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    plot!(URange, freqsOsc[mode,:], c=modeColors[mode], lw=lw,  label=false)
end
plt32 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.15,0.05])
for mode in 1:nModes
    plot!(URange, dampsOsc[mode,:]./freqsOsc[mode,:], c=modeColors[mode], lw=lw,  label="Mode $mode")
end
plt3 = plot(plt31,plt32, layout=(2,1))
display(plt3)
savefig(string(absPath,"/SMWFlutter_Vgf.pdf"))

# Mode shapes at specified velocity
plt4 = plot(xlabel="\$x_1/L\$", ylabel="\$u_2\$", title="Chordwise bending mode shape at U=$(URange[ind]) m/s")
plt5 = plot(xlabel="\$x_1/L\$", ylabel="\$u_3\$", title="Flapwise bending mode shape at U=$(URange[ind]) m/s")
plt6 = plot(xlabel="\$x_1/L\$", ylabel="\$p_1\$", title="Torsional mode shape at U=$(URange[ind]) m/s")
for m in 1:nModes
    plot!(plt4,x1_0/L, u2_modeShapes[m], c=modeColors[m], lw=lw, label="Mode $m")
    plot!(plt5,x1_0/L, u3_modeShapes[m], c=modeColors[m], lw=lw, label="Mode $m")
    plot!(plt6,x1_0/L, p1_modeShapes[m], c=modeColors[m], lw=lw, label="Mode $m")
end
display(plt4)
display(plt5)
display(plt6)
savefig(string(absPath,"/SMWFlutter_modalDisp.pdf"))

println("Finished SMWFlutter.jl")