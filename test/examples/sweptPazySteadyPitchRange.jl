using AeroBeams, DelimitedFiles

# Sweep angle [rad]
Λ = 30*π/180

# Root pitch angle range
θRange = π/180*vcat(0,3,5,7)

# Airspeed range
URange = collect(1:1:100)

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-undef"

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(flatPlate)

# Flag for upright position
upright = false

# Gravity
g = 0

# System solver
σ0 = 0.5
σstep = 0.5
maxIter = 50
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep,maximumIterations=maxIter)

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(θRange),length(URange))
tipTwist = Array{Float64}(undef,length(θRange),length(URange))
tipAoA = Array{Float64}(undef,length(θRange),length(URange))
AoA = Array{Vector{Float64}}(undef,length(θRange),length(URange))
cn = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u1_of_x1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u2_of_x1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))

# Sweep pitch angle
for (i,θ) in enumerate(θRange)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(θ*180/π)) deg, U = $U m/s")
        # Model
        model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,Λ=Λ,θ=θ,airspeed=U,g=g,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepStructuralCorrections=sweepStructuralCorrections)
        # Create and solve problem
        problem[i,j] = create_SteadyProblem(model=model,systemSolver=NR)
        solve!(problem[i,j])
        # Outputs
        tipOOP[i,j] = problem[i,j].nodalStatesOverσ[end][nElem].u_n2_b[3]
        tip_p = problem[i,j].nodalStatesOverσ[end][nElem].p_n2_b
        R = first(rotation_tensor_WM(tip_p))
        Δ = R*AeroBeams.a2
        tipTwist[i,j] = asind(Δ[3])
        tipAoA[i,j] = problem[i,j].model.elements[end].aero.flowAnglesAndRates.αₑ*180/π
        cn[i,j] = [problem[i,j].model.elements[e].aero.aeroCoefficients.cn for e in 1:nElem]
        AoA[i,j] = [problem[i,j].model.elements[e].aero.flowAnglesAndRates.αₑ for e in 1:nElem]
        u1_of_x1[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[1],problem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
        u2_of_x1[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[2],problem[i,j].nodalStatesOverσ[end][e].u_n2[2]) for e in 1:nElem]...)
        u3_of_x1[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[3],problem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
    end
end

# Undeformed nodal and elemental positions
x1_0 = vcat([vcat(problem[1].model.elements[e].r_n1[1],problem[1].model.elements[e].r_n2[1]) for e in 1:nElem]...)
x2_0 = vcat([vcat(problem[1].model.elements[e].r_n1[2],problem[1].model.elements[e].r_n2[2]) for e in 1:nElem]...)
x3_0 = vcat([vcat(problem[1].model.elements[e].r_n1[3],problem[1].model.elements[e].r_n2[3]) for e in 1:nElem]...)
x1_e = getindex.(getfield.(problem[1].model.elements, :x1), 1)

# Deformed nodal positions
x1_def = [x1_0 .+ u1_of_x1[i,j] for i in eachindex(θRange), j in eachindex(URange)]
x2_def = [x2_0 .+ u2_of_x1[i,j] for i in eachindex(θRange), j in eachindex(URange)]
x3_def = [x3_0 .+ u3_of_x1[i,j] for i in eachindex(θRange), j in eachindex(URange)]

# Load reference data (from AePW4 meetings)
dispΛ30θ5U60_Nastran = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/dispLambda30AoA5U60_Nastran.txt")

println("Finished sweptPazySteadyPitchRange.jl")