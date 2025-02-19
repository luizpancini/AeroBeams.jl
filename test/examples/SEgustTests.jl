using AeroBeams, DelimitedFiles

# Core function
function SEgustTestsCore(aeroSolver,gustLoadsSolver,testCase)

    # Set test case data
    if testCase == 1
        Ma = 0.2
        U = 68.06
        b = 0.40663
        w = 1.1878
        t₀ = 600*b/U
        tf = t₀ + 60*b/U
        θ = 0*π/180 
    elseif testCase == 2
        Ma = 0.2
        U = 68.06
        b = 0.40663
        w = 1.1878
        t₀ = 600*b/U
        tf = t₀ + 60*b/U
        θ = 10*π/180
    elseif testCase == 3
        Ma = 0.2
        U = 68.06
        b = 0.40663
        w = 1.1878
        t₀ = 600*b/U
        tf = t₀ + 60*b/U
        θ = 15*π/180   
    end

    # Gust
    gust = create_SharpEdgedGust(initialTime=t₀,verticalVelocity=w)

    # Wing surface
    airfoil = deepcopy(NACA0012)
    update_Airfoil_params!(airfoil,Ma=Ma,U=U,b=b)
    derivationMethod = AD()
    surf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=2*b,normSparPos=1/4,updateAirfoilParameters=true)

    # Wing beam
    L = 1
    nElem = 1
    ∞ = 1e12
    wing = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],rotationParametrization="E321",p0=[0;0;θ],aeroSurface=surf)

    # BCs
    clamp1 = create_BC(name="clamp1",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
    clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

    # Model
    SEgustTests = create_Model(name="SEgustTests",beams=[wing],BCs=[clamp1,clamp2],v_A=[0;U;0],gust=gust)

    # Set system solver options
    σ0 = 1.0
    maxIter = 100
    NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

    # Time variables
    Δt = (tf-t₀)/1000

    # Initial velocities update options
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=false, relaxFactor=0.5, Δt=Δt/10)

    # Create and solve dynamic problem
    problem = create_DynamicProblem(model=SEgustTests,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions,displayProgress=true)
    solve!(problem)

    # Unpack numerical solution
    t = problem.timeVector
    α = [problem.aeroVariablesOverTime[i][1].flowAnglesAndRates.αₑ for i in 1:length(t)]
    cn = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:length(t)]
    ct = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.ct for i in 1:length(t)]
    χ = [problem.elementalStatesOverTime[i][1].χ for i in 1:length(t)]
    if typeof(aeroSolver) == BLi
        fN = [problem.aeroVariablesOverTime[i][1].BLiFlow.fN for i in 1:length(t)]
        fPrimeN = [problem.aeroVariablesOverTime[i][1].BLiFlow.fPrimeN for i in 1:length(t)]
        f2PrimeN = [problem.elementalStatesOverTime[i][1].χ[4] for i in 1:length(t)]
    end
    cl = @. cn*cos(θ) + ct*sin(θ)

    # Load reference data
    if testCase == 1
        ΔclRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "gustTests", "SE_A0.txt"))
    elseif testCase == 2
        ΔclRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "gustTests", "SE_A10.txt"))
    elseif testCase == 3
        ΔclRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "gustTests", "SE_A15.txt"))
    end

    # Time index of gust encounter
    ind = argmin(abs.(t .- t₀))

    # Non-dimensional time and cl increment vectors
    τ = U/b * (t[ind:end] .- t[ind])
    Δcl = cl[ind:end] .- cl[ind]

    return τ,Δcl,ΔclRef
end

# Test cases
tests = 1:3

# Aerodynamic and gust indicial solvers
aeroSolvers = [QuasiSteady(); Indicial(); Inflow(); BLi(); BLo()]
gustLoadsSolvers = [IndicialGust("Kussner"); IndicialGust("Berci&Righi")]

# Initialize outputs
τ = Array{Vector{Float64}}(undef,length(aeroSolvers),length(gustLoadsSolvers),3)
Δcl = Array{Vector{Float64}}(undef,length(aeroSolvers),length(gustLoadsSolvers),3)
ΔclRef = Array{Matrix{Float64}}(undef,length(aeroSolvers),length(gustLoadsSolvers),3)

# Loop aerodynamic solver
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Loop gust solver
    for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
        # Loop test cases
        for (k,testCase) in enumerate(tests)
            # Display status
            println("Solver $i, gust solver $j, test $k")
            # Solve for current configuration
            τ[i,j,k],Δcl[i,j,k],ΔclRef[i,j,k] = SEgustTestsCore(aeroSolver,gustLoadsSolver,testCase)
        end
    end
end

println("Finished SEgustTests.jl")
