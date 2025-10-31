using AeroBeams, DelimitedFiles

# Hinge configurations
hingeConfigurations = ["free"; "locked"]

# Gust frequency range [Hz]
ωRange = unique(vcat(0.2:0.2:1,1:0.25:2,2:0.5:4,4:1:7,7:2:15))

# Gust maximum vertical velocity [m/s]
Ug = 0.5

# Root pitch angle [rad]
θ = 7.5*π/180

# Airspeed [m/s]
U = 18

# Stiffness of the spring around the hinge for in-plane bending
kIPBendingHinge = 1e4

# Gravity
g = 9.80665

# Discretization
nElementsInner = 16
nElementsFFWT = 4
nElem = nElementsInner + nElementsFFWT

# Tip loss options
hasTipCorrection = true
tipLossDecayFactor = 12

# Solution method for hinge constraint
solutionMethod = "addedResidual"

# Time variables
Δt = [2e-3, 5e-3]
tf = max.(2, 2 ./ ωRange)

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Initialize outputs
t = Array{Vector{Float64}}(undef,2,length(ωRange))
M2root = Array{Vector{Float64}}(undef,2,length(ωRange))
ΔM2peaks = Array{Vector{Float64}}(undef,2,length(ωRange))
problem = Array{DynamicProblem}(undef,2,length(ωRange))

# Loop hinge configurations
for (c,config) in enumerate(hingeConfigurations)
    # Loop gust frequency range
    for (i,ω) in enumerate(ωRange)
        # Display progress
        display("Solving for config = $(config), ω = $ω Hz")
        # Gust
        gust = create_OneMinusCosineGust(initialTime=0,duration=1/ω,verticalVelocity=Ug)
        # Model
        model = create_HealyBaselineFFWT(solutionMethod=solutionMethod,hingeConfiguration=config,kIPBendingHinge=kIPBendingHinge,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,g=g,gust=gust) 
        # Solve steady problem for initial conditions
        steadyProblem = create_SteadyProblem(model=model,systemSolver=NR)
        solve!(steadyProblem)
        # Create and solve dynamic problem
        problem[c,i] = create_DynamicProblem(model=model,finalTime=tf[i],Δt=Δt[c],systemSolver=NR,skipInitialStatesUpdate=true,x0=steadyProblem.x)
        solve!(problem[c,i])
        # Unpack numerical solution
        t[c,i] = problem[c,i].savedTimeVector
        M2root[c,i] = [problem[c,i].nodalStatesOverTime[j][1].M_n1[2] for j in 1:length(t[c,i])]
        ΔM2root = M2root[c,i] .- M2root[c,i][1]
        # Find peaks
        ΔM2peaksTime = [t[c,i][argmin(ΔM2root)]; t[c,i][argmax(ΔM2root)]]
        ΔM2peaks[c,i] = -[minimum(ΔM2root); maximum(ΔM2root)]
        println("Minimum and maximum ΔM2: $(ΔM2peaks[c,i])")
        # Plot
        plt = plot(xlabel="Time [s]", ylabel="ΔWRBM [N.m]", title="\$\\omega_g\$ = $ω Hz")
        plot!(t[c,i], -(M2root[c,i] .- M2root[c,i][1]), c=:black, lw=2, ls=:solid, label=false)
        scatter!([ΔM2peaksTime[1]], [ΔM2peaks[c,i][1]], c=:red, ms=3, msw=0, label=false)
        scatter!([ΔM2peaksTime[2]], [ΔM2peaks[c,i][2]], c=:red, ms=3, msw=0, label=false)
        display(plt)
    end
end

# Load reference data
ΔWRBM_locked_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTlockedVfreeOMCGustFloating/DWRBM_locked.txt")
ΔWRBM_free_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTlockedVfreeOMCGustFloating/DWRBM_free.txt")

println("Finished HealyBaselineFFWTlockedVfreeOMCGustFloating.jl")