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

# Stiffness factor for torsion
λT = 5

# Altitude
h = 20e3

# Discretization
nElemWing = 80
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# System solvers
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)
NReigen = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=0.5,displayStatus=false)

# Set number of vibration modes
nModes = 8

# Set bending curvature and airspeed ranges
k2Range = range(-0.015,0.045,5)
URange = unique(sort(vcat(20:0.5:60)))

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(k2Range),length(URange))
eigenProblem = Array{EigenProblem}(undef,length(k2Range),length(URange))

trimAoA = fill(NaN, length(k2Range), length(URange))

untrackedFreqs = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
untrackedDamps = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
untrackedEigenvectors = [fill(NaN64+im*NaN64, nModes, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
freqs = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
damps = [fill(NaN64, nModes) for k2 in 1:length(k2Range), U in 1:length(URange)]
modeDampings = [fill(NaN64, length(URange)) for k2 in 1:length(k2Range), mode in 1:nModes]
modeFrequencies = [fill(NaN64, length(URange)) for k2 in 1:length(k2Range), mode in 1:nModes]

highestConvUindex = Array{Int64}(undef,length(k2Range))

x1_0 = Array{Vector{Float64}}(undef, length(k2Range))
x3_0 = Array{Vector{Float64}}(undef, length(k2Range))
x1_n = Array{Vector{Float64}}(undef, length(k2Range))
x1_e = Array{Vector{Float64}}(undef, length(k2Range))
u1_of_x1 = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
x1_def = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
x3_def = Array{Vector{Float64}}(undef, length(k2Range),length(URange))

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for k2 = $k2, U = $U m/s")
        # Model for trim problem
        cHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,torsionStiffnessFactor=λT,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag,altitude=h)
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
        println("Trim AoA = $(trimAoA[i,j]*180/π)")
        # Model for eigen problem
        wingModel,_ = create_SMW(aeroSolver=aeroSolver,stiffnessFactor=λ,torsionStiffnessFactor=λT,airspeed=U,nElem=div(nElemWing,2),altitude=h,cd0=wingCd0,k2=k2,hasInducedDrag=hasInducedDrag,θ=trimAoA[i,j])
        # Set initial guess solution as previous known solution
        x0Eig = j == 1 ? zeros(0) : eigenProblem[i,j-1].x
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=wingModel,nModes=nModes,frequencyFilterLimits=[1,Inf],systemSolver=NReigen,x0=x0Eig)
        solve!(eigenProblem[i,j])
        # Skip if unconverged
        if !eigenProblem[i,j].systemSolver.convergedFinalSolution
            highestConvUindex[i] = j-1
            break
        else
            highestConvUindex[i] = j
        end
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = eigenProblem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(eigenProblem[i,j].dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = eigenProblem[i,j].eigenvectorsOscillatoryCplx
        # Undeformed jig-shape properties
        if j == 1
            # Undeformed nodal positions of right wing
            x1_0[i] = vcat([vcat(wingModel.elements[e].r_n1[1],wingModel.elements[e].r_n2[1]) for e in 1:div(nElemWing,2)]...)
            x3_0[i] = vcat([vcat(wingModel.elements[e].r_n1[3],wingModel.elements[e].r_n2[3]) for e in 1:div(nElemWing,2)]...)
            # Nodal and elemental arclength positions
            x1_n[i] = vcat([vcat(wingModel.elements[e].x1_n1,wingModel.elements[e].x1_n2) for e in 1:div(nElemWing,2)]...)
            x1_e[i] = [wingModel.elements[e].x1 for e in 1:div(nElemWing,2)]
        end
        # Displacements over span
        u1_of_x1[i,j] = vcat([vcat(eigenProblem[i,j].nodalStatesOverσ[end][e].u_n1[1],eigenProblem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:div(nElemWing,2)]...)
        u3_of_x1[i,j] = vcat([vcat(eigenProblem[i,j].nodalStatesOverσ[end][e].u_n1[3],eigenProblem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:div(nElemWing,2)]...)
        u1_of_x1[i,j] .-= u1_of_x1[i,j][1]
        u3_of_x1[i,j] .-= u3_of_x1[i,j][1]
        # Deformed nodal positions
        x1_def[i,j] = x1_0[i] .+ u1_of_x1[i,j]
        x3_def[i,j] = x3_0[i] .+ u3_of_x1[i,j]
    end
    # Frequencies and dampings after mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking_hungarian(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies and dampings by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URange)]
    end
end

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_flutter_matchedPitch_k2_range_torsionStiff.jl"
absPath = string(pwd(),relPath)
mkpath(absPath)

using Plots, ColorSchemes

# Plot configurations
colors = cgrad(:rainbow, length(k2Range), categorical=true)
ts = 10
fs = 16
lfs = 10
tsz = 10
lw = 2
ms = 3
msw = 0
L = 16
mshape = [:circle, :star, :utriangle, :pentagon, :diamond]
labels = ["\$k_2 = $(k2)\$" for k2 in k2Range]
lambdaTstr = string("lambdaT",λT)
gr()

# Root locus
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,2], ylims=[0,100], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
display(plt_RL)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_torsionStiff_",lambdaTstr,"_rootlocus.pdf"))

# Root locus - focus on 1st T-IP mode
plt_RL_TIP1 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-1,2], ylims=[0,30], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    for mode in 1:nModes
        plot!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        plot!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
display(plt_RL_TIP1)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_torsionStiff_",lambdaTstr,"_rootlocus_TIP1.pdf"))

# Root locus - focus on 2nd T-IP mode
plt_RL_TIP2 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-2,1], ylims=[50,90], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    for mode in 1:nModes
        plot!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        plot!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
display(plt_RL_TIP2)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_torsionStiff_",lambdaTstr,"_rootlocus_TIP2.pdf"))

# V-g-f
for (i,k2) in enumerate(k2Range)
    plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,60], tickfont=font(ts), guidefont=font(12))
    for mode in 1:nModes
        plot!(URange, modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
    plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.25,0.15], tickfont=font(ts), guidefont=font(12), legendfontsize=lfs, legend=:topleft)
    for mode in 1:nModes
        plot!(URange, modeDampings[i,mode]./modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
    plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
    display(plt_Vgf)
    savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_torsionStiff_",lambdaTstr,"_Vgf",i,".pdf"))
end

# Normalized deformed span at airspeed limits
plt_disp = plot(xlabel="Normalized spanwise direction", ylabel="Normalized vertical direction", xlims=[0,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
plot!([NaN], [NaN], ls=:solid, c=:black, lw=lw, label=string("\$U_{\\infty} = ",URange[1]," \$ m/s"))
plot!([NaN], [NaN], ls=:dashdot, c=:black, lw=lw, label=string("\$U_{\\infty} = ",URange[end]," \$ m/s"))
for (i,k2) in enumerate(k2Range)
    plot!(x1_def[i,1]/L, x3_def[i,1]/L, ls=:solid, c=colors[i], lw=lw, label=labels[i])
    plot!(x1_def[i,end]/L, x3_def[i,end]/L, ls=:dashdot, c=colors[i], lw=lw, label=false)
end
display(plt_disp)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_torsionStiff_",lambdaTstr,"_disp.pdf"))

println("Finished cHALEwing_flutter_matchedPitch_k2_range_torsionStiff.jl")