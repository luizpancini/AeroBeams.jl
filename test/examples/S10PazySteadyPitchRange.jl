using AeroBeams, DelimitedFiles

# Sweep angle [rad]
Λ = 10*π/180

# Root pitch angle offset (assumed)
θoffset = 0.0*π/180

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
tipMass = 11e-3
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
        tipOOP[i,j] = problem[i,j].nodalStatesOverσ[end][nElem].u_n2_b[3]
        tipAoA[i,j] = problem[i,j].model.elements[end].aero.flowAnglesAndRates.αₑ*180/π
        ϵ11Root = problem[i,j].compElementalStatesOverσ[end][1].γ[1]
        κ2Root = problem[i,j].compElementalStatesOverσ[end][1].κ[2]
        κ3Root = problem[i,j].compElementalStatesOverσ[end][1].κ[3]
        rootEps[i,j] = (ϵ11Root - κ2Root*zSG - κ3Root*ySG)
    end
end

# Load reference data (experiments from Technion - https://doi.org/10.5281/zenodo.16912763)
steady_strains_aoa0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_strains_aoa0.txt")
steady_strains_aoa1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_strains_aoa1.txt")
steady_strains_aoa3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_strains_aoa3.txt")
steady_strains_aoa5_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_strains_aoa5.txt")
steady_strains_aoa7_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_strains_aoa7.txt")
steady_strains_aoa10_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/steady_strains_aoa10.txt")

println("Finished S10PazySteadyPitchRange.jl")