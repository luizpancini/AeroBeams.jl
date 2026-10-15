using AeroBeams, JLD2, Interpolations, Plots, ColorSchemes

# Reference data IDs
refIDs = [5209]

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Flag to update tip correction with airspeed
tipLossFunctionIsAirspeedDependent = true

# Flag to smooth airspeed signal and corresponding moving average window
smoothAirspeedSignal = true
smoothingWindow = Int(1e3)

# Aerodynamic solver
aeroSolver = BLi()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Time variables
Δt = 2.5e-4
trackingFrequency = max(1,round(Int,2.5e-4/Δt))

# Set system solver options
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(maximumIterations=maxIter,relativeTolerance=relTol,displayStatus=false,alwaysUpdateJacobian=true,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2,allowAdvanceThroughUnconvergedAeroStates=aeroSolver.name=="BLi")

# Root strain gauge coordinates on the cross-section (the spar cs is 60 x 2.25 mm)
ySG_LE = 50e-3/2    # y-position: on LE
ySG_TE = -50e-3/2   # y-position: on TE
zSG = 2.25e-3/2     # z-position: on top

# Initialize outputs
dynamicProblem = Array{DynamicProblem}(undef,length(refIDs))
Uprofile = Array{Function}(undef,length(refIDs))
t = Array{Vector{Float64}}(undef,length(refIDs))
tipAoA = Array{Vector{Float64}}(undef,length(refIDs))
tipOOP = Array{Vector{Float64}}(undef,length(refIDs))
rootEpsLE = Array{Vector{Float64}}(undef,length(refIDs))
rootEpsTE = Array{Vector{Float64}}(undef,length(refIDs))

# Sweep IDs
for (i,refID) in enumerate(refIDs)
    # Load reference data (from Technion)
    @load pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/ID_"*string(refID)*"_flutter_data.jld2" data
    Λ = data["sweepAngle"]
    θ = data["pitchAngle"]
    tipMassConfig = data["tipMassConfig"]
    tipMass = data["tipMass"]
    tipMassOffset = data["tipMassOffset"]
    t_exp = vcat(data["t_ref"]...)
    U_exp = vcat(data["U_ref"]...)
    root_eps_LE_exp = vcat(data["root_eps_LE_ref"]...)
    root_eps_TE_exp = vcat(data["root_eps_TE_ref"]...)
    # tipMassOffset = 2e-2
    # println("Setting tipMassOffset to $(tipMassOffset)")
    # Display configuration
    println("Configuration: Λ=$(Λ*180/π) deg, θ=$(θ*180/π) deg, tip mass of $(tipMass*1e3) g, offset by $(tipMassOffset*1e2) cm of the $(tipMassConfig)")
    println("Runtime of $(round(t_exp[end]-t_exp[1],digits=1)) s, airspeed ranging from $(round(minimum(U_exp),digits=1)) to $(round(maximum(U_exp),digits=1)) m/s")
    # Set tip mass position
    tipMassPos = tipMassConfig == "LE" ? chord*normSparPos + tipMassOffset : -(chord*(1-normSparPos) + tipMassOffset)
    # Smooth airspeed signal, if applicable
    U_exp = smoothAirspeedSignal ? moving_average(U_exp,smoothingWindow) : U_exp
    # Airspeed profile (interpolate linearly between experimental time steps)
    itp = interpolate((t_exp,), U_exp, Gridded(Linear()))
    etp = extrapolate(itp, Flat())
    Uprofile[i] = t -> etp(t)
    # Model for steady problem
    steadyModel = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,Λ=Λ,θ=θ,airspeed=U_exp[1],hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepStructuralCorrections=sweepStructuralCorrections,tipMass=tipMass,ηtipMass=[0;tipMassPos;0]))
    # Create and solve steady problem for initial airspeed
    steadyProblem = create_SteadyProblem(model=steadyModel)
    solve!(steadyProblem)
    # Model for dynamic problem
    dynamicModel = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,Λ=Λ,θ=θ,airspeed=Uprofile[i],hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepStructuralCorrections=sweepStructuralCorrections,tipMass=tipMass,ηtipMass=[0;tipMassPos;0],tipLossFunctionIsAirspeedDependent=tipLossFunctionIsAirspeedDependent))
    # Create and solve dynamic problem
    dynamicProblem[i] = create_DynamicProblem(model=dynamicModel,initialTime=t_exp[1],finalTime=t_exp[end],Δt=Δt,systemSolver=NR,x0=steadyProblem.x,trackingFrequency=trackingFrequency)
    solve!(dynamicProblem[i])
    # Unpack numerical solution
    t[i] = dynamicProblem[i].savedTimeVector
    tipAoA[i] = [dynamicProblem[i].aeroVariablesOverTime[j][nElem].flowAnglesAndRates.αₑ for j in eachindex(t[i])]
    tipOOP[i] = -[dynamicProblem[i].nodalStatesOverTime[j][nElem].u_n2[1] for j in eachindex(t[i])]
    ϵ11Root = [dynamicProblem[i].compElementalStatesOverTime[j][1].γ[1] for j in eachindex(t[i])]
    κ2Root = [dynamicProblem[i].compElementalStatesOverTime[j][1].κ[2] for j in eachindex(t[i])]
    κ3Root = [dynamicProblem[i].compElementalStatesOverTime[j][1].κ[3] for j in eachindex(t[i])]
    rootEpsLE[i] = (ϵ11Root .- κ2Root*zSG .- κ3Root*ySG_LE)
    rootEpsTE[i] = (ϵ11Root .- κ2Root*zSG .- κ3Root*ySG_TE)
end

# Set paths
relPath = "/dev/sweptPazy/S10/outputs/S10PazyMatchingAirspeedSweep"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
showTitle = true
showAirspeedInPlots = false
interactivePlots = true
colors = cgrad(:rainbow, 2, categorical=true)
epsColors = cgrad(:rainbow, 3, categorical=true)
ts = 10
fs = 16
lfs = 8
lw = 2
ms = 6
msw = 1
α = 0.5
lsExp = :solid
if interactivePlots
    plotlyjs()
else
    gr()
end

# Plot limits
oopYLIMS = Dict(4023 => [ 0, 50],
                4024 => [25, 50],
                4032 => [ 0, 25],
                4033 => [ 0, 25],
                4034 => [15, 35],
                4035 => [15, 35],
                4303 => [20, 27.5],
                4306 => [20, 30],
                4307 => [20, 30],
                5208 => [15, 30],
                5209 => [20, 40])
aoaYLIMS = Dict(4023 => [-5, 25],
                4024 => [-5, 25],
                4032 => [-10, 25],
                4033 => [-10, 25],
                4034 => [-10, 25],
                4035 => [-10, 25],
                4303 => [ 2, 12],
                4306 => [ 0, 15],
                4307 => [ 0, 15],
                5208 => [-5, 25],
                5209 => [-5, 25])
epsYLIMS = Dict(4023 => [1500, 3000],
                4024 => [2000, 3200],
                4032 => [   0, 3200],
                4033 => [   0, 3200],
                4034 => [   0, 3200],
                4035 => [   0, 3200],
                4303 => [1000, 2300],
                4306 => [1000, 2400],
                4307 => [1000, 2400],
                5208 => [   0, 3200],
                5209 => [   0, 3200])

# Window for spectrogram
rootEpsWindow = Dict(4023 => 2^13,
                     4024 => 2^13,
                     4032 => 2^13,
                     4033 => 2^13,
                     4034 => 2^13,
                     4035 => 2^13,
                     4303 => 2^13,
                     4306 => 2^13,
                     4307 => 2^13,
                     5208 => 2^13,
                     5209 => 2^14)   

# Sweep IDs
for (i,refID) in enumerate(refIDs)
    # Case ID
    idString = string("ID", refID, "_", tipLossType, "_", aeroSolver.name)
    # Load reference data (from Technion)
    @load pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/ID_"*string(refID)*"_flutter_data.jld2" data
    Λ = data["sweepAngle"]
    θ = data["pitchAngle"]
    tipMassConfig = data["tipMassConfig"]
    tipMass = data["tipMass"]
    tipMassOffset = data["tipMassOffset"]
    t_exp = vcat(data["t_ref"]...)
    U_exp = vcat(data["U_ref"]...)
    root_eps_LE_exp = vcat(data["root_eps_LE_ref"]...)
    root_eps_TE_exp = vcat(data["root_eps_TE_ref"]...)
    # Title
    title = showTitle ? "Λ=$(Λ*180/π), θ=$(θ*180/π), TM of $(tipMass*1e3) g, offset $(tipMassOffset*1e2) cm of the $(tipMassConfig)" : ""
    # Airspeed
    plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]", title=title, xlims=extrema(t_exp), tickfont=font(ts), guidefont=font(fs))
    if smoothAirspeedSignal
        plot!(t_exp, U_exp, lw=lw, c=:blue, alpha=α, label=false)
    end
    plot!(t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    display(plt_U)
    savefig(string(absPath,"/S10PazyMatchingAirspeedSweep_",idString,"_U.pdf"))
    # Tip OOP displacement
    plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", title=title, xlims=extrema(t_exp), ylims=oopYLIMS[refID], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
    plot!(t[i], tipOOP[i]/L*100, lw=lw, c=colors[1], label=false)
    plot!([NaN], [NaN], lw=lw, ls=lsExp, c=colors[2], label=false)
    if showAirspeedInPlots
        plt_tipOOP_twin = twinx(plt_tipOOP)
        plot!(plt_tipOOP_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_tipOOP_twin, t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    end
    display(plt_tipOOP)
    savefig(string(absPath,"/S10PazyMatchingAirspeedSweep_",idString,"_tipOOP.pdf"))
    # Tip AoA
    plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", title=title, xlims=extrema(t_exp), ylims=aoaYLIMS[refID], tickfont=font(ts), guidefont=font(fs))
    plot!(t[i], tipAoA[i]*180/π, lw=lw, c=colors[1], label=false)
    plot!([NaN], [NaN], lw=lw, ls=lsExp, c=colors[2], label=false)
    if showAirspeedInPlots
        plt_tipAOA_twin = twinx(plt_tipAOA)
        plot!(plt_tipAOA_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_tipAOA_twin, t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    end
    display(plt_tipAOA)
    savefig(string(absPath,"/S10PazyMatchingAirspeedSweep_",idString,"_tipAoA.pdf"))
    # Root LE axial strains
    plt_rootEpsLE = plot(xlabel="Time [s]", ylabel="Root LE strains (\$\\mu\$)", title=title, xlims=extrema(t_exp), ylims=epsYLIMS[refID], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
    plot!(t[i], rootEpsLE[i]*1e6, lw=lw, c=epsColors[1], alpha=α, label="AeroBeams")
    plot!(t_exp, root_eps_LE_exp, lw=lw, ls=lsExp, c=epsColors[2], alpha=α, label="Revivo & Raveh (2025)")
    if showAirspeedInPlots
        plt_rootEpsLE_twin = twinx(plt_rootEpsLE)
        plot!(plt_rootEpsLE_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_rootEpsLE_twin, t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    end
    display(plt_rootEpsLE)
    savefig(string(absPath,"/S10PazyMatchingAirspeedSweep_",idString,"_rootEpsLE.pdf"))
    # Root TE axial strains
    plt_rootEpsTE = plot(xlabel="Time [s]", ylabel="Root TE strains (\$\\mu\$)", title=title, xlims=extrema(t_exp), ylims=epsYLIMS[refID], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
    plot!(t[i], rootEpsTE[i]*1e6, lw=lw, c=epsColors[1], alpha=α, label=false)
    plot!(t_exp, root_eps_TE_exp, lw=lw, ls=lsExp, c=epsColors[2], alpha=α, label=false)
    if showAirspeedInPlots
        plt_rootEpsTE_twin = twinx(plt_rootEpsTE)
        plot!(plt_rootEpsTE_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_rootEpsTE_twin, t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    end
    display(plt_rootEpsTE)
    savefig(string(absPath,"/S10PazyMatchingAirspeedSweep_",idString,"_rootEpsTE.pdf"))
    # Average root strains spectrogram
    rootEps = (rootEpsLE[i].+rootEpsTE[i])/2
    f_eps, t_eps, S_eps = spectrogram(t[i], rootEps*1e6, window=rootEpsWindow[refID])
    plt_eps_spectrum = heatmap(t_eps, f_eps, S_eps; ylims=[0,100], title=title, xlabel="Time [s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
    display(plt_eps_spectrum)
    savefig(string(absPath,"/S10PazyMatchingAirspeedSweep_",idString,"_rootEpsSpectrogram.pdf"))
    # Experimental average root strains spectrogram
    rootEps_exp = (root_eps_LE_exp.+root_eps_TE_exp)/2
    rootEps_exp = AeroBeams.fix_nans(rootEps_exp)
    f_eps, t_eps, S_eps = spectrogram(t_exp, rootEps_exp, window=2^12)
    plt_eps_spectrum = heatmap(t_eps, f_eps, S_eps; ylims=[0,100], title=title, xlabel="Time [s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
    display(plt_eps_spectrum)
    savefig(string(absPath,"/S10PazyMatchingAirspeedSweep_",idString,"_rootEpsSpectrogram_exp.pdf"))
end

println("Finished S10PazyMatchingAirspeedSweep.jl")