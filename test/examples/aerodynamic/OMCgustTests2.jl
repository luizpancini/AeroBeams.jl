using AeroBeams

# Core function
function OMCgustTestsCore(aeroSolver,gustSolver,testCase)

    # Set test case data
    if testCase == 1
        Ma = 0.2            # Mach number
        U = 68.06           # Airspeed [m/s]
        H = π               # Normalized gust length
        b = 0.40663         # Airfoil semichord [m]
        τ = 2*H*b/U         # Gust duration [s]
        w = 2.3748          # Gust peak velocity [m/s]
        t₀ = 60*b/U         # Time of gust encounter [s]
        tf = t₀ + 30*b/U    # Total simulation time [s]
        θ = 0*π/180         # Airfoil pitch angle [deg]
    elseif testCase == 2
        Ma = 0.2
        U = 68.06
        H = π
        b = 0.40663
        τ = 2*H*b/U
        w = 2.3748
        t₀ = 60*b/U
        tf = t₀ + 30*b/U
        θ = 10*π/180  
    elseif testCase == 3
        Ma = 0.2
        U = 68.06
        H = π
        b = 0.40663
        τ = 2*H*b/U
        w = 2.3748
        t₀ = 60*b/U
        tf = t₀ + 30*b/U
        θ = 15*π/180  
    elseif testCase == 4
        Ma = 0.2
        U = 68.06
        H = 8π
        b = 0.40663
        τ = 2*H*b/U
        w = 2.3748
        t₀ = 60*b/U
        tf = t₀ + 80*b/U
        θ = 0*π/180
    elseif testCase == 5
        Ma = 0.2
        U = 68.06
        H = 8π
        b = 0.40663
        τ = 2*H*b/U
        w = 2.3748
        t₀ = 60*b/U
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
    update_all_airfoil_parameters!(airfoil,Ma)
    derivationMethod = AD()
    surf = create_AeroSurface(ppSolver=aeroSolver,gustSolver=gustSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=2*b,normSparPos=1/4,updateAirfoilParameters=true)

    # Wing beam
    L = 1
    nElem = 1
    wing = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e12)],I=[inertia_matrix(ρA=1)],p0=transform_to_WM([0,0,θ],from=:E321),surface=surf)

    # BCs
    clamp1 = create_BC(name="clamp1",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
    clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

    # Model
    OMCgustTests = create_Model(name="OMCgustTests",beams=[wing],BCs=[clamp1,clamp2],v_A=[0,U,0],gust=gust)

    # Time variables
    timeSteps = 1500
    Δt = (tf-t₀)/timeSteps

    # Create and solve dynamic problem
    problem = create_DynamicProblem(model=OMCgustTests,finalTime=tf,Δt=Δt)
    solve!(problem)

    # Unpack numerical solution
    sol = problem.solutionHistory
    aeroVars = sol.aeroVariables
    t = sol.t
    cn = [aeroVars[i][1].aeroCoefficients.cn for i in eachindex(t)]
    ct = [aeroVars[i][1].aeroCoefficients.ct for i in eachindex(t)]
    cl = @. cn*cos(θ) + ct*sin(θ)

    # Load reference data
    if testCase == 1
        ΔclRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/gustTests/Hpi_A0.txt")
    elseif testCase == 2
        ΔclRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/gustTests/Hpi_A10.txt")
    elseif testCase == 3
        ΔclRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/gustTests/Hpi_A15.txt")
    elseif testCase == 4
        ΔclRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/gustTests/H8pi_A0.txt")
    elseif testCase == 5
        ΔclRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/gustTests/H8pi_A10.txt")
    elseif testCase == 6
        ΔclRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/gustTests/H8pi_A15.txt")
    end

    # Time index of gust encounter
    ind = argmin(abs.(t .- t₀)) - 1

    # Non-dimensional time and cl increment vectors
    τ = U/b * (t[ind:end] .- t[ind])
    Δcl = cl[ind:end] .- cl[ind]

    return τ,Δcl,ΔclRef

end 

# Tests range
tests = [3,6]

# Aerodynamic and gust indicial solvers
aeroSolvers = [QuasiSteady(); Indicial(); Inflow(); BLi(); BLo()]
gustLoadsSolvers = [IndicialGust(:Kussner)]

# Initialize outputs
τ = Array{Vector{Float64}}(undef,length(aeroSolvers),length(gustLoadsSolvers),length(tests))
Δcl = Array{Vector{Float64}}(undef,length(aeroSolvers),length(gustLoadsSolvers),length(tests))
ΔclRef = Array{Matrix{Float64}}(undef,length(aeroSolvers),length(gustLoadsSolvers),length(tests))

# Loop aerodynamic solver
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Loop gust solver
    for (j,gustSolver) in enumerate(gustLoadsSolvers)
        # Loop test cases
        for (k,testCase) in enumerate(tests)
            # Display status #src
            println("Solver $i, gust solver $j, test $testCase")
            τ[i,j,k],Δcl[i,j,k],ΔclRef[i,j,k] = OMCgustTestsCore(aeroSolver,gustSolver,testCase)
        end
    end
end

println("Finished OMCgustTests2.jl")
