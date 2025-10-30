using AeroBeams, LinearInterpolations

# Aerodynamic solver
aeroSolver = Indicial()

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Stiffness factor
λ = 1

# Altitude
h = 20e3

# Discretization
nElemWing = 80
nElemTailBoom = 5
nElemHorzStabilizer = 4
nElemVertStabilizer = 2

# System solvers
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,pseudoInverseMethod=:dampedLeastSquares,displayStatus=true)
NReigen = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=0.5,displayStatus=false)

# Set number of vibration modes
nModes = 8

# Linearly spanwise-varying bending curvature
k2root = 0.045
k2 = x1 -> k2root*(1-(x1/16)^(1))

# Set bending curvature and airspeed ranges
URange = vcat(20:0.5:50)

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(URange))
eigenProblem = Array{EigenProblem}(undef,length(URange))

trimAoA = fill(NaN, length(URange))

untrackedFreqs = [fill(NaN64, nModes) for U in 1:length(URange)]
untrackedDamps = [fill(NaN64, nModes) for U in 1:length(URange)]
untrackedEigenvectors = [fill(NaN64+im*NaN64, nModes, nModes) for U in 1:length(URange)]

x1_def = Array{Vector{Float64}}(undef,length(URange))
x3_def = Array{Vector{Float64}}(undef,length(URange))

# Sweep airspeed
for (j,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Model for trim problem
    cHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag,altitude=h)
    # Set initial guess solution as previous known solution
    x0Trim = j == 1 ? zeros(0) : trimProblem[j-1].x
    # Create and trim problem
    trimProblem[j] = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim,x0=x0Trim)
    solve!(trimProblem[j])
    # Extract trim variables
    trimAoA[j] = trimProblem[j].aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
    println("Trim AoA = $(trimAoA[j]*180/π)")
    # Model for eigen problem
    wingModel,_ = create_SMW(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElem=div(nElemWing,2),altitude=h,cd0=wingCd0,k2=k2,hasInducedDrag=hasInducedDrag,θ=trimAoA[j])
    # Set initial guess solution as previous known solution
    x0Eig = j == 1 ? zeros(0) : eigenProblem[j-1].x
    # Create and solve eigen problem
    eigenProblem[j] = create_EigenProblem(model=wingModel,nModes=nModes,frequencyFilterLimits=[1,Inf],systemSolver=NReigen,x0=x0Eig)
    solve!(eigenProblem[j])
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[j] = eigenProblem[j].frequenciesOscillatory
    untrackedDamps[j] = round_off!(eigenProblem[j].dampingsOscillatory,1e-8)
    untrackedEigenvectors[j] = eigenProblem[j].eigenvectorsOscillatoryCplx
    # Undeformed nodal positions of right wing
    x1_0 = vcat([vcat(wingModel.elements[e].r_n1[1],wingModel.elements[e].r_n2[1]) for e in 1:div(nElemWing,2)]...)
    x3_0 = vcat([vcat(wingModel.elements[e].r_n1[3],wingModel.elements[e].r_n2[3]) for e in 1:div(nElemWing,2)]...)
    # Displacements over span
    u1_of_x1 = vcat([vcat(eigenProblem[j].nodalStatesOverσ[end][e].u_n1[1],eigenProblem[j].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:div(nElemWing,2)]...)
    u3_of_x1 = vcat([vcat(eigenProblem[j].nodalStatesOverσ[end][e].u_n1[3],eigenProblem[j].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:div(nElemWing,2)]...)
    u1_of_x1 .-= u1_of_x1[1]
    u3_of_x1 .-= u3_of_x1[1]
    # Deformed nodal positions
    x1_def[j] = x1_0 .+ u1_of_x1
    x3_def[j] = x3_0 .+ u3_of_x1
end

# Frequencies and dampings after mode tracking
freqs,damps,_ = mode_tracking_hungarian(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

# Separate frequencies and dampings by mode
modeDampings = [fill(NaN64, length(URange)) for mode in 1:nModes]
modeFrequencies = [fill(NaN64, length(URange)) for mode in 1:nModes]
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[j][mode] for j in eachindex(URange)]
    modeDampings[mode] = [damps[j][mode] for j in eachindex(URange)]
end

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_flutter_matchedPitch_k2fun.jl"
absPath = string(pwd(),relPath)
mkpath(absPath)

using Plots, ColorSchemes

# Plot configurations
colors = palette([:royalblue, :blueviolet, :deeppink, :darkorange, :gold])
ts = 10
fs = 16
lfs = 10
tsz = 10
lw = 2
ms = 3
msw = 0
L = 16
k2str = string("k2",k2root)
gr()

# Root locus
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,2], ylims=[0,120], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=colors[mode], shape=:circle, ms=ms, msw=msw, label=false)
    scatter!([modeDampings[mode][1]], [modeFrequencies[mode][1]], c=colors[mode], shape=:circle, ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
end
display(plt_RL)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2fun_rootlocus_",k2str,".pdf"))

# V-g-f
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,60], tickfont=font(ts), guidefont=font(12))
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode], c=colors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.25,0.15], tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=:topleft)
for mode in 1:nModes
    plot!(URange, modeDampings[mode]./modeFrequencies[mode], c=colors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2fun_Vgf_",k2str,".pdf"))

# Normalized deformed span at airspeed limits
plt_disp = plot(xlabel="Normalized spanwise direction", ylabel="Normalized vertical direction", xlims=[0,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
plot!(x1_def[1]/L, x3_def[1]/L, ls=:solid, c=:black, lw=lw, label=string("\$U_{\\infty} = ",URange[1]," \$ m/s"))
plot!(x1_def[end]/L, x3_def[end]/L, ls=:dashdot, c=:black, lw=lw, label=string("\$U_{\\infty} = ",URange[end]," \$ m/s"))
display(plt_disp)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2fun_disp_",k2str,".pdf"))

println("Finished cHALEwing_flutter_matchedPitch_k2fun.jl")