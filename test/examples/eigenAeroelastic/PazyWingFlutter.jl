using AeroBeams

# Airspeed range [m/s]
URange = vcat(1e-2,0.5:0.5:100)

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Root pitch angle [rad]
θ = deg2rad(5)

# Initialize model
PazyWingFlutter,nElem,L,chord,normSparPos = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,pitchAngle=θ)

# Undeformed nodal and midpoint positions
wing = PazyWingFlutter.beams[1]
x1_0 = vcat([vcat(wing.elements[e].assemblyData.r_n1[3],wing.elements[e].assemblyData.r_n2[3]) for e in 1:nElem]...)
x3_0 = -vcat([vcat(wing.elements[e].assemblyData.r_n1[1],wing.elements[e].assemblyData.r_n2[1]) for e in 1:nElem]...)

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :arnoldi
nev = aeroSolver isa DynamicStallSolver ? 150 : 50

# Initialize outputs
problem = Array{EigenProblem}(undef,length(URange))
eigensols = Array{EigenSolution}(undef,length(URange))
tip_OOP = Array{Float64}(undef,length(URange))
u1 = Array{Vector{Float64}}(undef,length(URange))
u3 = Array{Vector{Float64}}(undef,length(URange))
x1_def = Array{Vector{Float64}}(undef,length(URange))
x3_def = Array{Vector{Float64}}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Update airspeed
    update_v_A!(PazyWingFlutter,v_A=[0,U,0])
    # Set initial guess solution
    x0 = i==1 ? zeros(0) : problem[i-1].systemArrays.x
    # Create and solve problem
    problem[i] = create_EigenProblem(model=PazyWingFlutter,eigenmethod=eigenmethod,nev=nev,x0=x0)
    solve!(problem[i])    
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
    # Tip OOP displacement at midchord
    states = problem[i].loadstepHistory.states
    tipTwist = twist_angle(states[end][nElem].nodal.p2_b, outputUnit=:deg)
    tip_OOP[i] = -(states[end][nElem].nodal.u2[1] - chord*(1/2-normSparPos)*sind(tipTwist))/L*100
    # Displacements over span
    u1[i] = vcat([vcat(states[end][e].nodal.u1_b[1],states[end][e].nodal.u2_b[1]) for e in 1:nElem]...)
    u3[i] = vcat([vcat(states[end][e].nodal.u1_b[3],states[end][e].nodal.u2_b[3]) for e in 1:nElem]...)
    # Deformed nodal positions
    x1_def[i] = x1_0 .+ u1[i]
    x3_def[i] = x3_0 .+ u3[i]
end

# Apply mode tracking
freqs,damps,modeShapes,_ = mode_tracking(eigensols)
rads2Hz!(freqs)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,modeShapesOsc,_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=1e-1, fmax=150, maxTracks=nModes)

# Flutter onset and offset speeds, frequencies, modes and tip OOP displacements
flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetModes,flutterOnsetDisp,flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetModes,flutterOffsetDisp = find_flutter(URange,freqsOsc,dampsOsc,tip_OOP)

# Show flutter results
for (m,speed,freq,disp) in zip(flutterOnsetModes,flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetDisp)
    println("Mode $m: onset speed = $(round(speed, digits=2)) m/s, frequency = $(round(freq, digits=2)) Hz, tip OOP = $(round(disp, digits=2)) %")
end
for (m,speed,freq,disp) in zip(flutterOffsetModes,flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetDisp)
    println("Mode $m: offset speed = $(round(speed, digits=2)) m/s, frequency = $(round(freq, digits=2)) Hz, tip OOP = $(round(disp, digits=2)) %")
end

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes
modesPlot = plot_mode_shapes(problem[end],modeShapes=modeShapesOsc[:,end],modeNumbers=collect(1:nModes),scale=0.2,legendPos=(0.25,0.98),view=(30,30),save=true,savePath=string(relPath,"/PazyWingFlutter_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 3
gr()

# Normalized deformed wingspan
plt0 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$", xlims=[0,1])
for (i,U) in enumerate(URange)
    plot!(x1_def[i]/L, x3_def[i]/L, lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt0)
savefig(string(absPath,"/PazyWingFlutter_disp.pdf"))

# V-g-f
plt11 = plot(ylabel="Frequency [Hz]")
for mode in 1:nModes
    plot!(URange, freqsOsc[mode,:], c=modeColors[mode], lw=lw,  label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.15,0.05], legend=:bottomleft)
for mode in 1:nModes
    plot!(URange, dampsOsc[mode,:]./Hz2rads.(freqsOsc[mode,:]), c=modeColors[mode], lw=lw,  label="Mode $mode")
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/PazyWingFlutter_Vgf.pdf"))

# Frequencies and dampings vs tip OOP displacement
plt21 = plot(ylabel="Frequency [Hz]")
for mode in 1:nModes
    plot!(tip_OOP, freqsOsc[mode,:], c=modeColors[mode], lw=lw,  label=false)
end
plt22 = plot(xlabel="Tip OOP displacement [% semispan]", ylabel="Damping Ratio", ylims=[-0.15,0.05], legend=:bottomleft)
for mode in 1:nModes
    plot!(tip_OOP, dampsOsc[mode,:]./Hz2rads.(freqsOsc[mode,:]), c=modeColors[mode], lw=lw,  label="Mode $mode")
end
plt2 = plot(plt21,plt22, layout=(2,1))
display(plt2)
savefig(string(absPath,"/PazyWingFlutter_OOPgf.pdf"))

println("Finished PazyWingFlutter.jl")