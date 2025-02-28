using AeroBeams, LinearInterpolations

# Mode tracking option
modeTracking = true

# Aerodynamic solver
aeroSolver = Indicial()

# Altitude
h = 20e3

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Discretization
nElemWing = 20
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# System solver for trim problem
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)

# Set number of vibration modes
nModes = 12

# Set bending curvature and airspeed ranges
k2Range = range(-0.015,0.045,5)
URange = collect(20:0.5:45)

# Initialize outputs
trimAoA = fill(NaN, length(k2Range), length(URange))
trimThrust = fill(NaN, length(k2Range), length(URange))
trimδ = fill(NaN, length(k2Range), length(URange))

untrackedFreqs = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
untrackedDamps = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
untrackedEigenvectors = [fill(NaN64+im*NaN64, nModes, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
freqs = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
damps = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
modeDampings = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
modeFrequencies = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
flutterOnsetSpeedOfMode = fill(NaN, length(k2Range), nModes)
flutterOffsetSpeedOfMode = fill(NaN, length(k2Range), nModes)
flutterOnsetSpeed = fill(NaN, length(k2Range))

x1_0 = Array{Vector{Float64}}(undef,length(k2Range))
x3_0 = Array{Vector{Float64}}(undef,length(k2Range))
x1_e_wing = Array{Vector{Float64}}(undef,length(k2Range))
u1_of_x1 = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
x1_def = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
x3_def = Array{Vector{Float64}}(undef,length(k2Range),length(URange))

highestConvUindex = Array{Int64}(undef,length(k2Range))

trimProblem = Array{TrimProblem}(undef,length(k2Range),length(URange))
eigenProblem = Array{EigenProblem}(undef,length(k2Range),length(URange))

# Attachment springs
μu = [1e-2; 1e-2; 1e-2; 1e-2; 1e-2]
μp = [10; 1e-2; 10; 1e-2; 1e-2]

# ELement ranges
elemRangeRightWing = 1 + div(nElemWing,2) : nElemWing

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Set attachment springs
    spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=μu[i]*[1; 1; 1],kp=μp[i]*[1; 1; 1])
    spring2 = create_Spring(elementsIDs=[nElemTailBoom],nodesSides=[2],ku=μu[i]*[1; 1; 1],kp=μp[i]*[1; 1; 1])
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for k2 = $k2, U = $U m/s")
        # Model for trim problem
        cHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
        # Add springs
        add_springs_to_beam!(beam=tailBoom,springs=[spring1,spring2])
        # Update model
        cHALEtrim.skipValidationMotionBasisA = true
        update_model!(cHALEtrim)
        # Set initial guess solution as previous known solution
        x0Trim = j == 1 ? zeros(0) : trimProblem[i,j-1].x
        # Create and trim problem
        trimProblem[i,j] = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim,x0=x0Trim)
        solve!(trimProblem[i,j])
        # Skip if unconverged
        if !trimProblem[i,j].systemSolver.convergedFinalSolution
            highestConvUindex[i] = j-1
            break
        else
            highestConvUindex[i] = j
        end
        # Extract trim variables
        trimAoA[i,j] = trimProblem[i,j].aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
        trimThrust[i,j] = stabilizersAero ? trimProblem[i,j].x[end-1]*trimProblem[i,j].model.forceScaling : trimProblem[i,j].x[end]*trimProblem[i,j].model.forceScaling
        trimδ[i,j] = stabilizersAero ? trimProblem[i,j].x[end] : 0
        println("Trim AoA = $(trimAoA[i,j]*180/π), trim thrust = $(trimThrust[i,j]), trim δ = $(trimδ[i,j]*180/π)")
        # Undeformed jig-shape properties
        if j == 1
            # Undeformed nodal positions of right wing
            x1_0[i] = vcat([vcat(cHALEtrim.elements[e].r_n1[1],cHALEtrim.elements[e].r_n2[1]) for e in elemRangeRightWing]...)
            x3_0[i] = vcat([vcat(cHALEtrim.elements[e].r_n1[3],cHALEtrim.elements[e].r_n2[3]) for e in elemRangeRightWing]...)
            # Undeformed elemental positions
            x1_e_wing[i] = [cHALEtrim.elements[e].x1 for e in elemRangeRightWing]
        end
        # Displacements over span
        u1_of_x1[i,j] = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].u_n1[1],trimProblem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in elemRangeRightWing]...)
        u3_of_x1[i,j] = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].u_n1[3],trimProblem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in elemRangeRightWing]...)
        u1_of_x1[i,j] .-= u1_of_x1[i,j][1]
        u3_of_x1[i,j] .-= u3_of_x1[i,j][1]
        # Deformed nodal positions
        x1_def[i,j] = x1_0[i] .+ u1_of_x1[i,j]
        x3_def[i,j] = x3_0[i] .+ u3_of_x1[i,j]
        # Model for eigen problem
        cHALEeigen,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElev=trimδ[i,j],thrust=trimThrust[i,j],k2=k2,hasInducedDrag=hasInducedDrag)
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=cHALEeigen,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],jacobian=trimProblem[i,j].jacobian[1:end,1:end-trimProblem[i,j].model.nTrimVariables],inertia=trimProblem[i,j].inertia,refTrimProblem=trimProblem[i,j])
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
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URange)]
    end
    # Flutter speed of each mode
    for mode in 1:nModes
        iOnset = findfirst(j -> modeDampings[i,mode][j] < 0 && modeDampings[i,mode][j+1] > 0, 1:length(URange)-1)
        iOffset = findfirst(j -> modeDampings[i,mode][j] > 0 && modeDampings[i,mode][j+1] < 0, 1:length(URange)-1)
        if isnothing(iOnset)
            flutterOnsetSpeedOfMode[i,mode] = Inf64
        else
            flutterOnsetSpeedOfMode[i,mode] = interpolate(modeDampings[i,mode][iOnset:iOnset+1],URange[iOnset:iOnset+1],0)
        end
        if isnothing(iOffset) || isnothing(iOnset)
            flutterOffsetSpeedOfMode[i,mode] = Inf64
        else
            flutterOffsetSpeedOfMode[i,mode] = interpolate(-modeDampings[i,mode][iOffset:iOffset+1],URange[iOffset:iOffset+1],0)
        end
    end
    flutterOnsetSpeed[i] = minimum(filter(!isinf,flutterOnsetSpeedOfMode[i,:]),init=Inf64)
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/outputs/figures/cHALE_flutter_k2_range"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at lowest airspeed
for (n,k2) in enumerate(k2Range)
    plt = plot_mode_shapes(eigenProblem[n,1],nModes=6,scale=5,view=(30,30),modalColorScheme=:rainbow,legendPos=:outertop,save=true,savePath=string(relPath,string("/cHALE_flutter_k2_range_modeShapesk2_",n,".pdf")))
    display(plt)
end

# Plot configurations
colors = cgrad(:rainbow, length(k2Range), categorical=true)
ts = 10
fs = 16
lfs = 10
tsz = 10
lw = 2
ms = 3
msw = 0
mshape = [:circle, :star, :utriangle, :pentagon, :diamond]
labels = ["\$k_2 = $(k2) \$" for k2 in k2Range]
L = 16
gr()

# Root locus
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-.03,0.02], ylims=[0.3,0.4], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:bottomleft)
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
display(plt_RL)
savefig(string(absPath,"/cHALE_flutter_k2_range_rootlocus.pdf"))

# Root locus (zoom on low frequency modes)
plt_RLlow = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,1.5], ylims=[0,20], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=(0.1,0.75))
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
phugoidTextPos = [0.75, 2.5]
TIPtextPos = [0.75, 7]
annotate!(phugoidTextPos[1], phugoidTextPos[2], text("phugoid", tsz))
annotate!(TIPtextPos[1], TIPtextPos[2], text("1st T-IP", tsz))
quiver!([phugoidTextPos[1],TIPtextPos[1]], [phugoidTextPos[2]-0.75,TIPtextPos[2]+0.5], quiver=([-0.5,-0.7], [-1.3,1.25]), arrow=:closed, linecolor=:black)
display(plt_RLlow)
savefig(string(absPath,"/cHALE_flutter_k2_range_rootlocus_low.pdf"))

# Root locus (zoom on T-IP mode)
plt_RLTIP = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-1.4,1.25], ylims=[3,21], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
display(plt_RLTIP)
savefig(string(absPath,"/cHALE_flutter_k2_range_rootlocus_TIP.pdf"))

# Root locus (zoom on phugoid mode)
plt_RLphugoid = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-0.15,0.15], ylims=[0,0.5], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
phugoidTextPos = [0.05, 0.45]
annotate!(phugoidTextPos[1], phugoidTextPos[2], text("phugoid", 15))
quiver!([phugoidTextPos[1]-0.03], [phugoidTextPos[2]], quiver=([-0.05], [-0.03]), arrow=:closed, linecolor=:black)
display(plt_RLphugoid)
savefig(string(absPath,"/cHALE_flutter_k2_range_rootlocus_phugoid.pdf"))

# V-g-f
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,50], tickfont=font(ts), guidefont=font(fs))
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.25,0.15], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    for mode in 1:nModes
        scatter!(plt_Vf, URange, modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        scatter!(plt_Vg, URange, modeDampings[i,mode]./modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,"/cHALE_flutter_k2_range_Vgf.pdf"))

# Flutter onset speeds vs k2
plt_Uf = plot(xlabel="\$k_2\$", ylabel="Flutter speed [m/s]", ylims=[0,45], xticks=k2Range, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:bottomleft)
plot!(k2Range,[22.42; 32.74; 22.42; 16.97; 12.50], marker=(:circle, 8), msw=0, c=:black, lw=lw, ls=:dash, label="Wing - zero pitch, no gravity")
plot!(k2Range,[41.49; 42.70; 39.96; 40.50; 20.0], marker=(:square, 8), msw=0, c=:black, lw=lw, ls=:dot, label="Wing - trim matched pitch")
plot!(k2Range,flutterOnsetSpeed, marker=(:diamond, 8), msw=0, c=:black, lw=lw, ls=:solid, label="Aircraft - trimmed flight")
display(plt_Uf)
savefig(string(absPath,"/cHALE_flutter_k2_range_speedOn.pdf"))

# Normalized deformed span at lowest and highest airspeeds
plt_u3 = plot(xlabel="Normalized spanwise direction", ylabel="Normalized OOP direction", xlims=[0,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], ls=:solid, c=:black, lw=lw, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
plot!([NaN], [NaN], ls=:dash, c=:black, lw=lw, label=string("\$U_{\\infty} = ",URange[end],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    plot!(x1_def[i,1]/L, x3_def[i,1]/L, ls=:solid, c=colors[i], lw=lw, label=false)
    plot!(x1_def[i,highestConvUindex[i]]/L, x3_def[i,highestConvUindex[i]]/L, ls=:dash, c=colors[i], lw=lw, label=false)
end
display(plt_u3)
savefig(string(absPath,"/cHALE_flutter_k2_range_disp.pdf"))

# Trim root angle of attack
plt_trimAoA = plot(xlabel="Airspeed [m/s]", ylabel="Wing root angle of attack [deg]", xlims=[URange[1],URange[end]], ylims=[-5,20], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimAoA[i,:]*180/π, c=colors[i], lw=lw, label=false)
end
display(plt_trimAoA)
savefig(string(absPath,"/cHALE_flutter_k2_range_trimAoA.pdf"))

# Trim thrust
plt_trimThrust = plot(xlabel="Airspeed [m/s]", ylabel="Thrust [N]", xlims=[URange[1],URange[end]], ylims=[0,70], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimThrust[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_trimThrust)
savefig(string(absPath,"/cHALE_flutter_k2_range_trimThrust.pdf"))

# Trim elevator deflection
plt_trimDelta = plot(xlabel="Airspeed [m/s]", ylabel="Elevator deflection [deg]", xlims=[URange[1],URange[end]], ylims=[-50,10], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimδ[i,:]*180/π, c=colors[i], lw=lw, label=false)
end
display(plt_trimDelta)
savefig(string(absPath,"/cHALE_flutter_k2_range_trimDelta.pdf"))

println("Finished cHALE_flutter_k2_range.jl")