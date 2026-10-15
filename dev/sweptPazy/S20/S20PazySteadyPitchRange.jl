using AeroBeams, DelimitedFiles, Plots, ColorSchemes

# Sweep angle [rad]
Λ = 20*π/180

# Root pitch angle offset (given by Technion)
θoffset = 0.4*π/180

# Root pitch angle range
θRange = π/180*vcat(0,1,3,5,7,10) .+ θoffset

# Airspeed range
URange = collect(1:1:70)

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-def"

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(flatPlate)

# Flag for upright position
upright = true

# Gravity
g = 9.80665

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Tip mass configuration (assumed at tipMassPosOffset behind TE)
tipMass = 15.5e-3
tipMassPosOffset = 0.05
tipMassPos = -(chord*(1-normSparPos) + tipMassPosOffset)

# Root strain gauge coordinates on the cross-section (the spar cs is 60 x 2.25 mm)
ySG = 0             # Averaged strains over LE and TE of spar
zSG = 2.25e-3/2     # On top of spar

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(θRange),length(URange))
tipAoA = Array{Float64}(undef,length(θRange),length(URange))
rootEps = Array{Float64}(undef,length(θRange),length(URange))

# Sweep pitch angle
for (i,θ) in enumerate(θRange)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(θ*180/π,digits=1)) deg, U = $U m/s")
        # Model
        model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,Λ=Λ,θ=θ,airspeed=U,g=g,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepStructuralCorrections=sweepStructuralCorrections,tipMass=tipMass,ηtipMass=[0;tipMassPos;0])
        # Create and solve problem
        problem[i,j] = create_SteadyProblem(model=model)
        solve!(problem[i,j])
        # Outputs
        tipOOP[i,j] = -problem[i,j].nodalStatesOverσ[end][nElem].u_n2[1]
        tipAoA[i,j] = problem[i,j].model.elements[end].aero.flowAnglesAndRates.αₑ*180/π
        ϵ11Root = problem[i,j].compElementalStatesOverσ[end][1].γ[1]
        κ2Root = problem[i,j].compElementalStatesOverσ[end][1].κ[2]
        κ3Root = problem[i,j].compElementalStatesOverσ[end][1].κ[3]
        rootEps[i,j] = (ϵ11Root - κ2Root*zSG - κ3Root*ySG)
    end
end

# Load reference data (experiments from Technion - https://doi.org/10.5281/zenodo.16354530)
steady_tipdisp_aoa0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipdisp_aoa0.txt")
steady_tipdisp_aoa1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipdisp_aoa1.txt")
steady_tipdisp_aoa3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipdisp_aoa3.txt")
steady_tipdisp_aoa5_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipdisp_aoa5.txt")
steady_tipdisp_aoa7_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipdisp_aoa7.txt")
steady_tipdisp_aoa10_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipdisp_aoa10.txt")

steady_tipaoa_aoa0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa0.txt")
steady_tipaoa_aoa1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa1.txt")
steady_tipaoa_aoa3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa3.txt")
steady_tipaoa_aoa5_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa5.txt")
steady_tipaoa_aoa7_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa7.txt")
steady_tipaoa_aoa10_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_tipaoa_aoa10.txt")

steady_strains_aoa0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_strains_aoa0.txt")
steady_strains_aoa1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_strains_aoa1.txt")
steady_strains_aoa3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_strains_aoa3.txt")
steady_strains_aoa5_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_strains_aoa5.txt")
steady_strains_aoa7_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_strains_aoa7.txt")
steady_strains_aoa10_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/steady_strains_aoa10.txt")

# Set paths
relPath = "/dev/sweptPazy/S20/outputs/S20PazySteadyPitchRange"
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
    plot!(URange, tipOOP[i,:]/L*100, c=colors[i], lw=lw, ls=:solid, label="\$\\theta = $(round(Int,θ*180/pi)) ^\\circ\$")
    if θ ≈ 0 + θoffset
        scatter!(steady_tipdisp_aoa0_ref[1,:], steady_tipdisp_aoa0_ref[2,:]*100, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 1*π/180 + θoffset
        scatter!(steady_tipdisp_aoa1_ref[1,:], steady_tipdisp_aoa1_ref[2,:]*100, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 3*π/180 + θoffset
        scatter!(steady_tipdisp_aoa3_ref[1,:], steady_tipdisp_aoa3_ref[2,:]*100, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 5*π/180 + θoffset
        scatter!(steady_tipdisp_aoa5_ref[1,:], steady_tipdisp_aoa5_ref[2,:]*100, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 7*π/180 + θoffset
        scatter!(steady_tipdisp_aoa7_ref[1,:], steady_tipdisp_aoa7_ref[2,:]*100, ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 10*π/180 + θoffset
        scatter!(steady_tipdisp_aoa10_ref[1,:], steady_tipdisp_aoa10_ref[2,:]*100, ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_tipOOP)
savefig(string(absPath,"/S20PazySteadyPitchRange_tipOOP.pdf"))

# Tip AoA vs. airspeed
plt_tipAOA = plot(xlabel="Airspeed [m/s]", ylabel="Tip angle of attack [deg]", xlims=[0,70], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tipAoA[i,:], c=colors[i], lw=lw, label=false)
    if θ ≈ 0 + θoffset
        scatter!(steady_tipaoa_aoa0_ref[1,:], steady_tipaoa_aoa0_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 1*π/180 + θoffset
        scatter!(steady_tipaoa_aoa1_ref[1,:], steady_tipaoa_aoa1_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 3*π/180 + θoffset
        scatter!(steady_tipaoa_aoa1_ref[1,:], steady_tipaoa_aoa3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 5*π/180 + θoffset
        scatter!(steady_tipaoa_aoa5_ref[1,:], steady_tipaoa_aoa5_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 7*π/180 + θoffset
        scatter!(steady_tipaoa_aoa7_ref[1,:], steady_tipaoa_aoa7_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 10*π/180 + θoffset
        scatter!(steady_tipaoa_aoa10_ref[1,:], steady_tipaoa_aoa10_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_tipAOA)
savefig(string(absPath,"/S20PazySteadyPitchRange_tipAoA.pdf"))

# Root axial strains vs. airspeed
plt_rootEps = plot(xlabel="Airspeed [m/s]", ylabel="Root axial strains (\$\\mu\$)", xlims=[0,70], tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
for (i,θ) in enumerate(θRange)
    plot!(URange, rootEps[i,:]*1e6, c=colors[i], lw=lw, ls=:solid, label=false)
    if θ ≈ 0 + θoffset
        scatter!(steady_strains_aoa0_ref[1,:], steady_strains_aoa0_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 1*π/180 + θoffset
        scatter!(steady_strains_aoa1_ref[1,:], steady_strains_aoa1_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 3*π/180 + θoffset
        scatter!(steady_strains_aoa3_ref[1,:], steady_strains_aoa3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 5*π/180 + θoffset
        scatter!(steady_strains_aoa5_ref[1,:], steady_strains_aoa5_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 7*π/180 + θoffset
        scatter!(steady_strains_aoa7_ref[1,:], steady_strains_aoa7_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif θ ≈ 10*π/180 + θoffset
        scatter!(steady_strains_aoa10_ref[1,:], steady_strains_aoa10_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_rootEps)
savefig(string(absPath,"/S20PazySteadyPitchRange_rootEps.pdf"))

println("Finished S20PazySteadyPitchRange.jl")