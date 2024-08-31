using AeroBeams, DelimitedFiles

# Core function
function OMCgustTestsCore(aeroSolver,gustLoadsSolver,testCase)

    # Set test case data
    if testCase == 1
        Ma = 0.2
        U = 68.06
        H = π
        b = 0.40663
        τ = 2*H*b/U
        w = 2.3748
        t₀ = 80*b/U
        tf = t₀ + 30*b/U
        θ = 0*π/180 
    elseif testCase == 2
        Ma = 0.2
        U = 68.06
        H = π
        b = 0.40663
        τ = 2*H*b/U
        w = 2.3748
        t₀ = 80*b/U
        tf = t₀ + 30*b/U
        θ = 10*π/180  
    elseif testCase == 3
        Ma = 0.2
        U = 68.06
        H = π
        b = 0.40663
        τ = 2*H*b/U
        w = 2.3748
        t₀ = 80*b/U
        tf = t₀ + 30*b/U
        θ = 15*π/180  
    elseif testCase == 4
        Ma = 0.2
        U = 68.06
        H = 8π
        b = 0.40663
        τ = 2*H*b/U
        w = 2.3748
        t₀ = 80*b/U
        tf = t₀ + 80*b/U
        θ = 0*π/180
    elseif testCase == 5
        Ma = 0.2
        U = 68.06
        H = 8π
        b = 0.40663
        τ = 2*H*b/U
        w = 2.3748
        t₀ = 80*b/U
        tf = t₀ + 80*b/U
        θ = 10*π/180
    elseif testCase == 6
        Ma = 0.2
        U = 68.06
        H = 8π
        b = 0.40663
        τ = 2*H*b/U
        w = 2.3748
        t₀ = 60*b/U
        tf = t₀ + 80*b/U
        θ = 15*π/180         
    end

    # Gust
    gust = create_OneMinusCosineGust(initialTime=t₀,duration=τ,verticalVelocity=w)

    # Wing surface
    airfoil = deepcopy(NACA0012)
    surf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=2*b,normSparPos=1/4,updateAirfoilParameters=true)

    # Wing beam
    L = 1
    nElem = 1
    ∞ = 1e12
    wing = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],rotationParametrization="E321",p0=[0;0;θ],aeroSurface=surf)

    # BCs
    clamp1 = create_BC(name="clamp1",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
    clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

    # Model
    OMCgustTests = create_Model(name="OMCgustTests",beams=[wing],BCs=[clamp1,clamp2],v_A=[0;U;0],gust=gust)

    # Set system solver options
    σ0 = 1.0
    maxIter = 100
    NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

    # Time variables
    Δt = (tf-t₀)/1000

    # Initial velocities update options
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=false, relaxFactor=0.5, Δt=Δt/10)

    # Create and solve dynamic problem
    problem = create_DynamicProblem(model=OMCgustTests,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
    solve!(problem)

    # Unpack numerical solution
    t = problem.timeVector
    α = [problem.aeroVariablesOverTime[i][1].flowAnglesAndRates.αₑ for i in 1:length(t)]
    cn = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:length(t)]
    ct = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.ct for i in 1:length(t)]
    if typeof(aeroSolver) == BLi
        fN = [problem.aeroVariablesOverTime[i][1].BLiFlow.fN for i in 1:length(t)]
        fPrimeN = [problem.aeroVariablesOverTime[i][1].BLiFlow.fPrimeN for i in 1:length(t)]
        f2PrimeN = [problem.elementalStatesOverTime[i][1].χ[4] for i in 1:length(t)]
    end
    cl = @. cn*cos(θ) + ct*sin(θ)

    # Load reference data
    if testCase == 1
        ΔclRef = readdlm("test/referenceData/gustTests/Hpi_A0.txt")
    elseif testCase == 2
        ΔclRef = readdlm("test/referenceData/gustTests/Hpi_A10.txt")
    elseif testCase == 3
        ΔclRef = readdlm("test/referenceData/gustTests/Hpi_A15.txt")
    elseif testCase == 4
        ΔclRef = readdlm("test/referenceData/gustTests/H8pi_A0.txt")
    elseif testCase == 5
        ΔclRef = readdlm("test/referenceData/gustTests/H8pi_A10.txt")
    elseif testCase == 6
        ΔclRef = readdlm("test/referenceData/gustTests/H8pi_A15.txt")
    end

    # Time index of gust encounter
    ind = argmin(abs.(t .- t₀))

    # Non-dimensional time and cl increment vectors
    τ = U/b * (t[ind:end] .- t[ind])
    Δcl = cl[ind:end] .- cl[ind]

    return τ,Δcl,ΔclRef
end      

# Aerodynamic and gust indicial solvers
aeroSolvers = [QuasiSteady(); Indicial(); Inflow(); BLi(); BLo()]
gustLoadsSolvers = [IndicialGust("Kussner"); IndicialGust("Berci&Righi")]

# Initialize outputs
τ = Array{Vector{Float64}}(undef,length(aeroSolvers),length(gustLoadsSolvers),6)
Δcl = Array{Vector{Float64}}(undef,length(aeroSolvers),length(gustLoadsSolvers),6)
ΔclRef = Array{Matrix{Float64}}(undef,length(aeroSolvers),length(gustLoadsSolvers),6)

# Loop aerodynamic solver
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Loop gust solver
    for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
        # Loop test cases
        for (k,testCase) in enumerate(1:6)
            # Display status
            println("Solver $i, gust solver $j, test $k")
            # Solve for current configuration
            τ[i,j,k],Δcl[i,j,k],ΔclRef[i,j,k] = OMCgustTestsCore(aeroSolver,gustLoadsSolver,testCase)
        end
    end
end

println("Finished OMCgustTests.jl")
