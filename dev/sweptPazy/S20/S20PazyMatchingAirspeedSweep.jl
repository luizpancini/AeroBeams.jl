using AeroBeams, JLD2, Interpolations, Plots, ColorSchemes

# Reference data IDs
refIDs = vcat(3020,3031,3034:3038,3201:3208)

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Flag to update tip correction with airspeed
tipLossFunctionIsAirspeedDependent = true

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Time variables
Δt = 5e-4
trackingFrequency = round(Int,1e-3/Δt)

# Set system solver options
maxIter = 50
relTol = 1e-9
NR = create_NewtonRaphson(maximumIterations=maxIter,relativeTolerance=relTol,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2,allowAdvanceThroughUnconvergedAeroStates=aeroSolver.name=="BLi")

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
    @load pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/ID_"*string(refID)*"_flutter_data.jld2" data
    Λ = data["sweepAngle"]
    θ = data["pitchAngle"]
    tipMassConfig = data["tipMassConfig"]
    tipMass = data["tipMass"]
    tipMassOffset = data["tipMassOffset"]
    t_exp = vcat(data["t_ref"]...)
    U_exp = vcat(data["U_ref"]...)
    root_eps_LE_exp = vcat(data["root_eps_LE_ref"]...)
    root_eps_TE_exp = vcat(data["root_eps_TE_ref"]...)
    # Display configuration
    println("Configuration: Λ=$(Λ*180/π) deg, θ=$(θ*180/π) deg, tip mass of $(tipMass*1e3) g, offset by $(tipMassOffset*1e2) cm of the $(tipMassConfig)")
    println("Runtime of $(round(t_exp[end]-t_exp[1],digits=1)) s, airspeed ranging from $(round(minimum(U_exp),digits=1)) to $(round(maximum(U_exp),digits=1)) m/s")
    # Set tip mass position
    tipMassPos = tipMassConfig == "LE" ? chord*normSparPos + tipMassOffset : -(chord*(1-normSparPos) + tipMassOffset)
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
    tipAoA[i] = [dynamicProblem[i].aeroVariablesOverTime[i][nElem].flowAnglesAndRates.αₑ for i in eachindex(t[i])]
    tipOOP[i] = -[dynamicProblem[i].nodalStatesOverTime[i][nElem].u_n2[1] for i in eachindex(t[i])]
    ϵ11Root = [dynamicProblem[i].compElementalStatesOverTime[i][1].γ[1] for i in eachindex(t[i])]
    κ2Root = [dynamicProblem[i].compElementalStatesOverTime[i][1].κ[2] for i in eachindex(t[i])]
    κ3Root = [dynamicProblem[i].compElementalStatesOverTime[i][1].κ[3] for i in eachindex(t)]
    rootEpsLE[i] = (ϵ11Root .- κ2Root*zSG .- κ3Root*ySG_LE)
    rootEpsTE[i] = (ϵ11Root .- κ2Root*zSG .- κ3Root*ySG_TE)
end

# Set paths
relPath = "/dev/sweptPazy/S20/outputs/S20PazyMatchingAirspeedSweep"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
showTitle = true
showAirspeedInPlots = true
interactivePlots = false
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
oopYLIMS = Dict(3020 => [ 0, 50],
                3031 => [25, 50],
                3034 => [ 0, 25],
                3035 => [ 0, 25],
                3036 => [15, 35],
                3037 => [15, 35],
                3038 => [20, 27.5],
                3201 => [20, 30],
                3202 => [20, 30],
                3203 => [15, 30],
                3204 => [15, 30],
                3205 => [20, 30],
                3206 => [20, 30],
                3207 => [15, 30],
                3208 => [15, 30])
aoaYLIMS = Dict(3020 => [-5, 25],
                3031 => [-5, 25],
                3034 => [-10, 25],
                3035 => [-10, 25],
                3036 => [-10, 25],
                3037 => [-10, 25],
                3038 => [ 2, 12],
                3201 => [ 0, 15],
                3202 => [ 0, 15],
                3203 => [-5, 25],
                3204 => [-5, 25],
                3205 => [ 0, 15],
                3206 => [ 0, 15],
                3207 => [-5, 25],
                3208 => [-5, 25])
epsYLIMS = Dict(3020 => [1500, 3000],
                3031 => [2000, 3200],
                3034 => [   0, 3200],
                3035 => [   0, 3200],
                3036 => [   0, 3200],
                3037 => [   0, 3200],
                3038 => [1000, 2300],
                3201 => [1000, 2400],
                3202 => [1000, 2400],
                3203 => [   0, 3200],
                3204 => [   0, 3200],
                3205 => [1000, 2400],
                3206 => [1000, 2400],
                3207 => [   0, 3200],
                3208 => [   0, 3200])

# Sweep IDs
for (i,refID) in enumerate(refIDs)
    # Case ID
    idString = string("ID", refID, "_", tipLossType, "_", aeroSolver.name)
    # Load reference data (from Technion)
    @load pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/ID_"*string(refID)*"_flutter_data.jld2" data
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
    title = showTitle ? "Configuration: Λ=$(Λ*180/π) deg, θ=$(θ*180/π) deg, tip mass of $(tipMass*1e3) g, offset by $(tipMassOffset*1e2) cm of the $(tipMassConfig)" : ""
    # Airspeed
    plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]", title=title, xlims=extrema(t_exp), tickfont=font(ts), guidefont=font(fs))
    plot!(t[i], Uprofile[i].(t[i]), lw=lw, c=:black, label=false)
    display(plt_U)
    savefig(string(absPath,"/S20PazyMatchingAirspeedSweep_",idString,"_U.pdf"))
    # Tip OOP displacement
    plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", title=title, xlims=extrema(t_exp), ylims=oopYLIMS[refID], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
    plot!(t[i], tipOOP[i]/L*100, lw=lw, c=colors[1], label=false)
    plot!([NaN], [NaN], lw=lw, ls=lsExp, c=colors[2], label=false)
    if showAirspeedInPlots
        plt_tipOOP_twin = twinx(plt_tipOOP)
        plot!(plt_tipOOP_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_tipOOP_twin, t_exp, U_exp, lw=lw, c=:black, label=false)
    end
    display(plt_tipOOP)
    savefig(string(absPath,"/S20PazyMatchingAirspeedSweep_",idString,"_tipOOP.pdf"))
    # Tip AoA
    plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", title=title, xlims=extrema(t_exp), ylims=aoaYLIMS[refID], tickfont=font(ts), guidefont=font(fs))
    plot!(t[i], tipAoA[i]*180/π, lw=lw, c=colors[1], label=false)
    plot!([NaN], [NaN], lw=lw, ls=lsExp, c=colors[2], label=false)
    if showAirspeedInPlots
        plt_tipAOA_twin = twinx(plt_tipAOA)
        plot!(plt_tipAOA_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_tipAOA_twin, t_exp, U_exp, lw=lw, c=:black, label=false)
    end
    display(plt_tipAOA)
    savefig(string(absPath,"/S20PazyMatchingAirspeedSweep_",idString,"_tipAoA.pdf"))
    # Root LE axial strains
    plt_rootEpsLE = plot(xlabel="Time [s]", ylabel="Root LE strains (\$\\mu\$)", title=title, xlims=extrema(t_exp), ylims=epsYLIMS[refID], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
    plot!(t[i], rootEpsLE[i]*1e6, lw=lw, c=epsColors[1], alpha=α, label="AeroBeams")
    plot!(t_exp, root_eps_LE_exp, lw=lw, ls=lsExp, c=epsColors[2], alpha=α, label="Revivo & Raveh (2025)")
    if showAirspeedInPlots
        plt_rootEpsLE_twin = twinx(plt_rootEpsLE)
        plot!(plt_rootEpsLE_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_rootEpsLE_twin, t_exp, U_exp, lw=lw, c=:black, label=false)
    end
    display(plt_rootEpsLE)
    savefig(string(absPath,"/S20PazyMatchingAirspeedSweep_",idString,"_rootEpsLE.pdf"))
    # Root TE axial strains
    plt_rootEpsTE = plot(xlabel="Time [s]", ylabel="Root TE strains (\$\\mu\$)", title=title, xlims=extrema(t_exp), ylims=epsYLIMS[refID], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
    plot!(t[i], rootEpsTE[i]*1e6, lw=lw, c=epsColors[1], alpha=α, label=false)
    plot!(t_exp, root_eps_TE_exp, lw=lw, ls=lsExp, c=epsColors[2], alpha=α, label=false)
    if showAirspeedInPlots
        plt_rootEpsTE_twin = twinx(plt_rootEpsTE)
        plot!(plt_rootEpsTE_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t_exp))
        plot!(plt_rootEpsTE_twin, t_exp, U_exp, lw=lw, c=:black, label=false)
    end
    display(plt_rootEpsTE)
    savefig(string(absPath,"/S20PazyMatchingAirspeedSweep_",idString,"_rootEpsTE.pdf"))
    # Average root strains spectrogram
    rootEps = (rootEpsLE[i].+rootEpsTE[i])/2
    f_eps, t_eps, S_eps = spectrogram(t[i], rootEps*1e6, window=2^12)
    U_vec = Uprofile[i].(t[i])
    U_eps = U_vec[searchsortedfirst.(Ref(t[i]), t_eps)]
    plt_eps_spectrum = heatmap(t_eps, f_eps, S_eps; ylims=[0,100], title=title, xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
    display(plt_eps_spectrum)
    savefig(string(absPath,"/S20PazyMatchingAirspeedSweep_",idString,"_rootEpsSpectrogram.pdf"))
end

println("Finished S20PazyMatchingAirspeedSweep.jl")