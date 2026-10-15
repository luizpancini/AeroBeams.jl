using AeroBeams

# Core function
function SEgustTestsCore(aeroSolver,gustSolver,testCase)

    # Set test case data
    if testCase == 1
        Ma = 0.2
        U = 68.06
        b = 0.40663
        w = 1.1878
        t₀ = 60*b/U
        tf = t₀ + 60*b/U
        θ = deg2rad(0)
    elseif testCase == 2
        Ma = 0.2
        U = 68.06
        b = 0.40663
        w = 1.1878
        t₀ = 60*b/U
        tf = t₀ + 60*b/U
        θ = deg2rad(10)
    elseif testCase == 3
        Ma = 0.2
        U = 68.06
        b = 0.40663
        w = 1.1878
        t₀ = 60*b/U
        tf = t₀ + 60*b/U
        θ = deg2rad(15)
    end

    # Gust
    gust = create_SharpEdgedGust(initialTime=t₀,verticalVelocity=w)

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
    SEgustTests = create_Model(name="SEgustTests",beams=[wing],BCs=[clamp1,clamp2],v_A=[0,U,0],gust=gust)

    # Time variables
    timeSteps = 1000
    Δt = (tf-t₀)/timeSteps

    # Create and solve dynamic problem
    problem = create_DynamicProblem(model=SEgustTests,finalTime=tf,Δt=Δt,displayProgress=true)
    solve!(problem)

    # Unpack numerical solution
    sol = problem.solutionHistory
    aeroVars = sol.aeroVariables
    t = sol.t
    cn = [aeroVars[i][1].aeroCoefficients.cn for i in eachindex(t)]
    ct = [aeroVars[i][1].aeroCoefficients.ct for i in eachindex(t)]
    cl = @. cn*cos(θ) + ct*sin(θ)

    # Load reference data
    ΔclRef = if testCase == 1
        readdlm(pkgdir(AeroBeams,"test","referenceData","gustTests", "SE_A0.txt"))
    elseif testCase == 2
        readdlm(pkgdir(AeroBeams,"test","referenceData","gustTests", "SE_A10.txt"))
    elseif testCase == 3
        readdlm(pkgdir(AeroBeams,"test","referenceData","gustTests", "SE_A15.txt"))
    end

    # Time index of gust encounter
    ind = argmin(abs.(t .- t₀)) - 1

    # Non-dimensional time and cl increment vectors
    τ = U/b * (t[ind:end] .- t[ind])
    Δcl = cl[ind:end] .- cl[ind]

    return τ,Δcl,ΔclRef

end

# Test cases
tests = 1:3

# Aerodynamic and gust indicial solvers
aeroSolvers = [QuasiSteady(); Indicial(); Inflow(); BLi(); BLo()]
gustLoadsSolvers = [IndicialGust(:Kussner); IndicialGust(:BerciRighi)]

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
            # Display status
            println("Solver $i, gust solver $j, test $k")
            # Solve for current configuration
            τ[i,j,k],Δcl[i,j,k],ΔclRef[i,j,k] = SEgustTestsCore(aeroSolver,gustSolver,testCase)
        end
    end
end

## Plots
# Plot configurations
lw = 2
ms = 5
gr()

# Loop aerodynamic solver
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Loop gust solver
    for (j,gustSolver) in enumerate(gustLoadsSolvers)
        # Loop test cases
        for (k,testCase) in enumerate(tests)
            # Set paths
            relPath = string("/test/outputs/figures/SEgustTests/SEgustTests_",aeroSolver.name,"_",string(gustSolver.indicialFunction),"_test",testCase)
            absPath = string(pwd(),relPath)
            mkpath(absPath)
            # Lift coefficient increment over time
            plt1 = plot(xlabel="\$\\tau\$ [semichords]", ylabel="\$\\Delta c_l\$")
            plot!(τ[i,j,k], Δcl[i,j,k], color=:black, lw=lw, label="AeroBeams")
            scatter!(ΔclRef[i,j,k][1,:], ΔclRef[i,j,k][2,:], color=:black, ms=ms, label="Mallik & Raveh (2019)")
            display(plt1)
            savefig(string(absPath,"/SEgustTests_dcl.pdf"))
        end
    end
end

println("Finished SEgustTests.jl")