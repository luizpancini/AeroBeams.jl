using AeroBeams, JLD2

# Sweep angle [rad]
Λ = deg2rad(10)

# Root pitch angle offset (assumed) [rad]
θoffset = deg2rad(0.0)

# Root pitch angle range [rad]
θRange = deg2rad.([0,1,3,5,7,10]) .+ θoffset

# Airspeed range [m/s]
URange = vcat(1e-3,1:1:70)

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(FlatPlate)

# Flag for upright position
upright = true

# Geometric properties
nElem,L,chord,normSparPos,nodalPositions = geometrical_properties_Pazy(Λ=Λ)

# Tip mass configuration
tipMassConfig = :TE

# Root strain gauge coordinates on the cross-section (the spar cs is 60 x 2.25 mm)
ySG = -50e-3/2      # TE strains
zSG = 2.23e-3/2     # On top of spar

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(θRange),length(URange))
tipAoA = Array{Float64}(undef,length(θRange),length(URange))
rootEps = Array{Float64}(undef,length(θRange),length(URange))

# Sweep pitch angle
for (i,θ) in enumerate(θRange)
    # Model
    model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,sweepAngle=Λ,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweptTipMassConfig=tipMassConfig)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(rad2deg(θ),digits=1)) deg, U = $U m/s")
        # Update airspeed
        update_v_A!(model, v_A=[0,U,0])
        # Create and solve problem
        problem[i,j] = create_SteadyProblem(model=model,x0=x0)
        solve!(problem[i,j])
        # Final states
        state = problem[i,j].loadstepHistory.states[end]
        aeroVars = problem[i,j].loadstepHistory.aeroVariables[end]
        # Outputs
        tipOOP[i,j] = -state[nElem].nodal.u2[1]
        tipAoA[i,j] = rad2deg(aeroVars[nElem].unsteadyKin.αₑ)
        ϵ11Root = state[1].midpoint.γ[1]
        κ2Root = state[1].midpoint.κ[2]
        κ3Root = state[1].midpoint.κ[3]
        rootEps[i,j] = (ϵ11Root - κ2Root*zSG - κ3Root*ySG)
        # Update initial guess for next airspeed
        x0 = problem[i,j].systemArrays.x
    end
end

# Load static results derived from raw data (https://doi.org/10.5281/zenodo.16912763)
s = load(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/static_results.jld2")
tipdisp_aoa0_ref = permutedims(hcat(s["V"][1:6], s["tip"][1:6,3]))
tipdisp_aoa1_ref = permutedims(hcat(s["V"][7:12], s["tip"][7:12,3]))
tipdisp_aoa3_ref = permutedims(hcat(s["V"][13:18], s["tip"][13:18,3]))
tipdisp_aoa5_ref = permutedims(hcat(s["V"][19:24], s["tip"][19:24,3]))
tipdisp_aoa7_ref = permutedims(hcat(s["V"][25:30], s["tip"][25:30,3]))
tipdisp_aoa10_ref = permutedims(hcat(s["V"][31:end], s["tip"][31:end,3]))

tipAoA_aoa0_ref = permutedims(hcat(s["V"][1:6], 0 .+ s["twist"][1:6]))
tipAoA_aoa1_ref = permutedims(hcat(s["V"][7:12], 1 .+ s["twist"][7:12]))
tipAoA_aoa3_ref = permutedims(hcat(s["V"][13:18], 3 .+ s["twist"][13:18]))
tipAoA_aoa5_ref = permutedims(hcat(s["V"][19:24], 5 .+ s["twist"][19:24]))
tipAoA_aoa7_ref = permutedims(hcat(s["V"][25:30], 7 .+ s["twist"][25:30]))
tipAoA_aoa10_ref = permutedims(hcat(s["V"][31:end], 10 .+ s["twist"][31:end]))

rootTEε_aoa0_ref = permutedims(hcat(s["V"][1:6], s["TE"][1:6]))
rootTEε_aoa1_ref = permutedims(hcat(s["V"][7:12], s["TE"][7:12]))
rootTEε_aoa3_ref = permutedims(hcat(s["V"][13:18], s["TE"][13:18]))
rootTEε_aoa5_ref = permutedims(hcat(s["V"][19:24], s["TE"][19:24]))
rootTEε_aoa7_ref = permutedims(hcat(s["V"][25:30], s["TE"][25:30]))
rootTEε_aoa10_ref = permutedims(hcat(s["V"][31:end], s["TE"][31:end]))

# Load tip angle of attack data (digitized from https://doi.org/10.2514/1.J066240, seems more accurate than derived data)
steady_tipaoa_aoa0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_tipaoa_aoa0.txt")
steady_tipaoa_aoa1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_tipaoa_aoa1.txt")
steady_tipaoa_aoa3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_tipaoa_aoa3.txt")
steady_tipaoa_aoa5_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_tipaoa_aoa5.txt")
steady_tipaoa_aoa7_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_tipaoa_aoa7.txt")
steady_tipaoa_aoa10_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_tipaoa_aoa10.txt")

## Plots
# Set paths
relPath = "/dev/Pazy/S10/outputs/S10PazySteadyPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 4
msw = 0
gr()

# Tip OOP displacement vs. airspeed
plt_tipOOP = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,70], tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
plot!([NaN], [NaN], lw=lw, c=:black, label="AeroBeams")
scatter!([NaN], [NaN], ms=ms, msw=msw, c=:black, label="Exp.")
for (i,θ) in enumerate(θRange)
    plot!(URange, tipOOP[i,:]/L*100, c=colors[i], lw=lw, ls=:solid, label="\$\\theta = $(round(Int,rad2deg(θ))) ^\\circ\$")
    if θ ≈ 0 + θoffset
        scatter!(tipdisp_aoa0_ref[1,:], tipdisp_aoa0_ref[2,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(1) + θoffset
        scatter!(tipdisp_aoa1_ref[1,:], tipdisp_aoa1_ref[2,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(3) + θoffset
        scatter!(tipdisp_aoa3_ref[1,:], tipdisp_aoa3_ref[2,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(5) + θoffset
        scatter!(tipdisp_aoa5_ref[1,:], tipdisp_aoa5_ref[2,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(7) + θoffset
        scatter!(tipdisp_aoa7_ref[1,:], tipdisp_aoa7_ref[2,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(10) + θoffset
        scatter!(tipdisp_aoa10_ref[1,:], tipdisp_aoa10_ref[2,:]*1e-3*100/L, ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_tipOOP)
savefig(string(absPath,"/S10PazySteadyPitchRange_tipOOP.pdf"))

# Tip AoA vs. airspeed
plt_tipAOA = plot(xlabel="Airspeed [m/s]", ylabel="Tip angle of attack [deg]", xlims=[0,70], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tipAoA[i,:], c=colors[i], lw=lw, label=false)
    if θ ≈ 0 + θoffset
        scatter!(steady_tipaoa_aoa0_ref[1,:], steady_tipaoa_aoa0_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(1) + θoffset
        scatter!(steady_tipaoa_aoa1_ref[1,:], steady_tipaoa_aoa1_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(3) + θoffset
        scatter!(steady_tipaoa_aoa1_ref[1,:], steady_tipaoa_aoa3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(5) + θoffset
        scatter!(steady_tipaoa_aoa5_ref[1,:], steady_tipaoa_aoa5_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(7) + θoffset
        scatter!(steady_tipaoa_aoa7_ref[1,:], steady_tipaoa_aoa7_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(10) + θoffset
        scatter!(steady_tipaoa_aoa10_ref[1,:], steady_tipaoa_aoa10_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_tipAOA)
savefig(string(absPath,"/S10PazySteadyPitchRange_tipAoA.pdf"))

# Root axial strains at the trailing edge vs. airspeed
plt_rootEps = plot(xlabel="Airspeed [m/s]", ylabel="Root axial strains (\$\\mu\$)", xlims=[0,70], tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
for (i,θ) in enumerate(θRange)
    plot!(URange, rootEps[i,:]*1e6, c=colors[i], lw=lw, ls=:solid, label=false)
    if θ ≈ deg2rad(0) + θoffset
        scatter!(rootTEε_aoa0_ref[1,:], rootTEε_aoa0_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(1) + θoffset
        scatter!(rootTEε_aoa1_ref[1,:], rootTEε_aoa1_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(3) + θoffset
        scatter!(rootTEε_aoa3_ref[1,:], rootTEε_aoa3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(5) + θoffset
        scatter!(rootTEε_aoa5_ref[1,:], rootTEε_aoa5_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(7) + θoffset
        scatter!(rootTEε_aoa7_ref[1,:], rootTEε_aoa7_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ deg2rad(10) + θoffset
        scatter!(rootTEε_aoa10_ref[1,:], rootTEε_aoa10_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_rootEps)
savefig(string(absPath,"/S10PazySteadyPitchRange_rootEps.pdf"))

println("Finished S10PazySteadyPitchRange.jl")