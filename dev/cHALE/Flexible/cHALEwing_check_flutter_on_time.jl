using AeroBeams, LinearInterpolations

# Aerodynamic solver
aeroSolver = Indicial()

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Altitude
h = 20e3

# Discretization
nElemWing = 40
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# Set number of vibration modes
nModes = 25

# Bending pre-curvature
k2 = 0.0

# Set airspeed range
URange = collect(20:1:60)

# Time simulation variables
Δt = 5e-2
tf = 180

# System solvers
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)
NReigen = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(URange))
eigenProblem = Array{EigenProblem}(undef,length(URange))

trimAoA = fill(NaN, length(URange))
trimThrust = fill(NaN, length(URange))
trimδ = fill(NaN, length(URange))

untrackedFreqs = [fill(NaN64, nModes) for U in 1:length(URange)]
untrackedDamps = [fill(NaN64, nModes) for U in 1:length(URange)]
untrackedEigenvectors = [fill(NaN64+im*NaN64, nModes, nModes) for U in 1:length(URange)]
freqs = [fill(NaN64, nModes) for U in 1:length(URange)]
damps = [fill(NaN64, nModes) for U in 1:length(URange)]
modeDampings = [fill(NaN64, length(URange)) for mode in 1:nModes]
modeFrequencies = [fill(NaN64, length(URange)) for mode in 1:nModes]
flutterOnsetSpeedOfMode = fill(NaN, nModes)
flutterOffsetSpeedOfMode = fill(NaN, nModes)

x1_def = Array{Vector{Float64}}(undef,length(URange))
x3_def = Array{Vector{Float64}}(undef,length(URange))

# Sweep airspeed
for (j,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Model for trim problem
    cHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag,altitude=h)
    # Set initial guess solution as previous known solution
    x0Trim = j == 1 ? zeros(0) : trimProblem[j-1].x
    # Create and trim problem
    trimProblem[j] = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim,x0=x0Trim)
    solve!(trimProblem[j])
    # Skip if unconverged
    if !trimProblem[j].systemSolver.convergedFinalSolution
        break
    end
    # Extract trim variables
    trimAoA[j] = trimProblem[j].aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
    trimThrust[j] = trimProblem[j].x[end-1]*trimProblem[j].model.forceScaling
    trimδ[j] = trimProblem[j].x[end]
    println("Trim AoA = $(trimAoA[j]*180/π)")
    # Model for eigen problem
    wingModel,_ = create_SMW(aeroSolver=aeroSolver,airspeed=U,nElem=nElemWing,altitude=h,cd0=wingCd0,k2=k2,hasInducedDrag=hasInducedDrag,θ=trimAoA[j])
    # Set initial guess solution as previous known solution
    x0Eig = j == 1 ? zeros(0) : eigenProblem[j-1].x
    # Create and solve eigen problem
    eigenProblem[j] = create_EigenProblem(model=wingModel,nModes=nModes,frequencyFilterLimits=[0,Inf64],systemSolver=NReigen,x0=x0Eig)
    solve!(eigenProblem[j])
    # Skip if unconverged
    if !eigenProblem[j].systemSolver.convergedFinalSolution
        break
    end
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[j] = eigenProblem[j].frequenciesOscillatory
    untrackedDamps[j] = round_off!(eigenProblem[j].dampingsOscillatory,1e-8)
    untrackedEigenvectors[j] = eigenProblem[j].eigenvectorsOscillatoryCplx
    # Undeformed nodal positions of right wing
    x1_0 = vcat([vcat(wingModel.elements[e].r_n1[1],wingModel.elements[e].r_n2[1]) for e in 1:nElemWing]...)
    x3_0 = vcat([vcat(wingModel.elements[e].r_n1[3],wingModel.elements[e].r_n2[3]) for e in 1:nElemWing]...)
    # Displacements over span
    u1_of_x1 = vcat([vcat(eigenProblem[j].nodalStatesOverσ[end][e].u_n1[1],eigenProblem[j].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElemWing]...)
    u3_of_x1 = vcat([vcat(eigenProblem[j].nodalStatesOverσ[end][e].u_n1[3],eigenProblem[j].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElemWing]...)
    u1_of_x1 .-= u1_of_x1[1]
    u3_of_x1 .-= u3_of_x1[1]
    # Deformed nodal positions
    x1_def[j] = x1_0 .+ u1_of_x1
    x3_def[j] = x3_0 .+ u3_of_x1
end

# Frequencies and dampings after mode tracking
freqs,damps,_ = mode_tracking_hungarian(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and dampings by mode
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[j][mode] for j in eachindex(URange)]
    modeDampings[mode] = [damps[j][mode] for j in eachindex(URange)]
end

# Flutter speed of each mode and flutter speed indices
indicesFlutterOnset = fill(0,0)
for mode in 1:nModes
    iOnset = findfirst(j -> modeDampings[mode][j] < 0 && modeDampings[mode][j+1] > 0, 1:length(URange)-1)
    if modeDampings[mode][1] > 0
        flutterOnsetSpeedOfMode[mode] = URange[1]
        push!(indicesFlutterOnset, 1)
    elseif isnothing(iOnset)
        flutterOnsetSpeedOfMode[mode] = Inf64
    else
        flutterOnsetSpeedOfMode[mode] = interpolate(modeDampings[mode][iOnset:iOnset+1],URange[iOnset:iOnset+1],0)
        push!(indicesFlutterOnset, iOnset)
    end
end

# Add indices post-flutter and sort
append!(indicesFlutterOnset, indicesFlutterOnset .+ 1)
sort!(unique!(indicesFlutterOnset))

# All flutter onset speeds, in order
flutterOnsetSpeedsAll = sort(filter(!isinf,flutterOnsetSpeedOfMode))
println("Flutter onset speeds: ", flutterOnsetSpeedsAll)

# Initialize dynamic problem outputs
dynamicProblem = Array{DynamicProblem}(undef,length(indicesFlutterOnset))
t = Array{Vector{Float64}}(undef,length(indicesFlutterOnset))
wingAoA = Array{Vector{Float64}}(undef,length(indicesFlutterOnset))
airspeed = Array{Vector{Float64}}(undef,length(indicesFlutterOnset))
u3tip = Array{Vector{Float64}}(undef,length(indicesFlutterOnset))

# Sweep flutter onset speeds
for (ind,j) in enumerate(indicesFlutterOnset)
    println("Running dynamic for U = $(URange[j]) m/s")
    # Set checked elevator deflection profile
    Δδ = -1*π/180
    tδinit = 0.5
    tδramp = 0.5
    tδpeak = tδinit+tδramp
    tδfinal = tδpeak+tδramp
    δ = t -> ifelse(
        t <= tδinit, 
        trimδ[j],
        ifelse(
            t <= tδpeak, 
            trimδ[j] + Δδ * ((t-tδinit) / (tδpeak-tδinit)),
            ifelse(
                t <= tδfinal, 
                trimδ[j] + Δδ - Δδ * ((t-tδpeak) / (tδfinal-tδpeak)),
                trimδ[j]
            )
        )
    )
    # Model for dynamic problem
    cHALEdynamic,_ = create_conventional_HALE(aeroSolver=aeroSolver,altitude=h,airspeed=URange[j],nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,hasInducedDrag=hasInducedDrag,k2=k2,δElev=δ,thrust=trimThrust[j])
    # Create and solve dynamic problem
    dynamicProblem[ind] = create_DynamicProblem(model=cHALEdynamic,x0=trimProblem[j].x[1:end-2],finalTime=tf,Δt=Δt,skipInitialStatesUpdate=true)
    solve!(dynamicProblem[ind])
    # Get wing root and tip elements
    lTipElem = 1
    lRootElem = div(nElemWing,2)
    rRootElem = lRootElem+1
    # Get outputs
    t[ind] = dynamicProblem[ind].savedTimeVector
    wingAoA[ind] = [(dynamicProblem[ind].aeroVariablesOverTime[i][lRootElem].flowAnglesAndRates.αₑ+dynamicProblem[ind].aeroVariablesOverTime[i][rRootElem].flowAnglesAndRates.αₑ)/2 for i in 1:length(t[ind])]
    airspeed[ind] = [(dynamicProblem[ind].aeroVariablesOverTime[i][lRootElem].flowVelocitiesAndRates.U∞+dynamicProblem[ind].aeroVariablesOverTime[i][rRootElem].flowVelocitiesAndRates.U∞)/2 for i in 1:length(t[ind])]
    u3tip[ind] = [dynamicProblem[ind].nodalStatesOverTime[i][lTipElem].u_n1[3]/16*100 for i in 1:length(t[ind])]
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_check_flutter_on_time.jl"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = palette([:royalblue, :blueviolet, :deeppink, :darkorange, :gold], nModes)
colors = palette([:royalblue, :blueviolet, :deeppink, :darkorange, :gold], length(indicesFlutterOnset))
ts = 10
fs = 16
lfs = 10
tsz = 10
lw = 2
ms = 3
msw = 0
gr()

# V-g-f
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,41], tickfont=font(ts), guidefont=font(12))
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.2,0.1], tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=:topleft)
for mode in 1:nModes
    plot!(URange, modeDampings[mode]./modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,"/cHALEwing_check_flutter_on_time_Vgf.pdf"))

# Root AoA
plt_AoA = plot(xlabel="Time [s]", ylabel="Normalized root AoA", tickfont=font(ts), guidefont=font(fs))
for (ind,j) in enumerate(indicesFlutterOnset)
    plot!(t[ind], (wingAoA[ind]/wingAoA[ind][1]), c=colors[ind], lw=lw, label=string("\$U = \$",URange[j]," m/s"))
end
display(plt_AoA)
savefig(string(absPath,"/cHALEwing_check_flutter_on_time_AoA.pdf"))

# Airspeed
plt_U = plot(xlabel="Time [s]", ylabel="Normalized airspeed", tickfont=font(ts), guidefont=font(fs))
for (ind,j) in enumerate(indicesFlutterOnset)
    plot!(t[ind], airspeed[ind]/airspeed[ind][1], c=colors[ind], lw=lw, label=string("\$U = \$",URange[j]," m/s"))
end
display(plt_U)
savefig(string(absPath,"/cHALEwing_check_flutter_on_time_U.pdf"))

# Tip OOP disp
plt_u3 = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", tickfont=font(ts), guidefont=font(fs))
for (ind,j) in enumerate(indicesFlutterOnset)
    plot!(t[ind], u3tip[ind], c=colors[ind], lw=lw, label=string("\$U = \$",URange[j]," m/s"))
end
display(plt_u3)
savefig(string(absPath,"/cHALEwing_check_flutter_on_time_u3tip.pdf"))


println("Finished cHALEwing_check_flutter_on_time.jl")