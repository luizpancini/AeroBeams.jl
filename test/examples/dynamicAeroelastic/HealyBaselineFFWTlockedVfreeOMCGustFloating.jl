using AeroBeams

# Hinge configurations
hingeConfigurations = [:free, :locked]

# Gust frequency range [Hz]
ωRange = unique(vcat(0.2:0.2:1,1:0.25:2,2:0.5:4,4:1:7,7:2:15))

# Gust maximum vertical velocity [m/s]
Ug = 0.5

# Root pitch angle [rad]
θ = deg2rad(7.5)

# Airspeed [m/s]
U = 18

# Discretization
nElementsInner = 16
nElementsFFWT = 4
nElem = nElementsInner + nElementsFFWT

# Tip loss options
hasTipCorrection = true
tipLossDecayFactor = 12

# Time variables
Δt = 5e-3
tf = max.(2, 2 ./ ωRange)

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Initialize outputs
t = Array{Vector{Float64}}(undef,2,length(ωRange))
M2root = Array{Vector{Float64}}(undef,2,length(ωRange))
ΔM2peaks = Array{Vector{Float64}}(undef,2,length(ωRange))

# Loop hinge configurations
for (c,config) in enumerate(hingeConfigurations)
    # Loop gust frequency range
    for (i,ω) in enumerate(ωRange)
        # Display progress
        display("Solving for config = $(config), ω = $ω Hz")
        # Gust
        gust = create_OneMinusCosineGust(duration=1/ω,verticalVelocity=Ug)
        # Model
        model = create_HealyBaselineFFWT(hingeConfiguration=config,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT,gust=gust) 
        # Solve steady problem for initial conditions
        steadyProblem = create_SteadyProblem(model=model,systemSolver=systemSolver)
        solve!(steadyProblem)
        # Create and solve dynamic problem
        dynProblem = create_DynamicProblem(model=model,finalTime=tf[i],Δt=Δt,systemSolver=systemSolver,x0=steadyProblem.systemArrays.x)
        solve!(dynProblem)
        # Unpack numerical solution
        sol = dynProblem.solutionHistory
        states = sol.states
        t[c,i] = sol.t
        M2root[c,i] = [states[j][1].nodal.M1[2] for j in eachindex(t[c,i])]
        ΔM2root = M2root[c,i] .- M2root[c,i][1]
        # Find peak values
        ΔM2peaks[c,i] = -[minimum(ΔM2root), maximum(ΔM2root)]
        # # Plot
        # plt = plot(xlabel="Time [s]", ylabel="ΔWRBM [N.m]", title="\$\\omega_g\$ = $ω Hz")
        # plot!(t[c,i], -(M2root[c,i] .- M2root[c,i][1]), c=:black, lw=2, ls=:solid, label=false)
        # scatter!([ΔM2peaksTime[1]], [ΔM2peaks[c,i][1]], c=:red, ms=3, msw=0, label=false)
        # scatter!([ΔM2peaksTime[2]], [ΔM2peaks[c,i][2]], c=:red, ms=3, msw=0, label=false)
        # display(plt)
    end
end

# Load reference data
ΔWRBM_locked_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTlockedVfreeOMCGustFloating/DWRBM_locked.txt")
ΔWRBM_free_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTlockedVfreeOMCGustFloating/DWRBM_free.txt")


## Plots
# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTlockedVfreeOMCGustFloating"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, 2, categorical=true)
ts = 12
fs = 16
lfs = 12
lw = 2
ms = 4
msw = 0
gr()

# Root OOP bending moment increment
plt_ΔWRBM = plot(xlabel="\$\\omega_g\$ [Hz]", ylabel="ΔWRBM [N.m]", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, xlims=[0,15], legend=(0.25,0.65))
plot!([NaN],[NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN],[NaN], c=:black, ms=ms, msw=msw, label="Healy (2023) - Nastran")
plot!([NaN],[NaN], c=colors[1], lw=lw, shape=:circle, ms=ms, msw=msw, label="Free")
plot!([NaN],[NaN], c=colors[2], lw=lw, shape=:circle, ms=ms, msw=msw, label="Locked")
for (c,config) in enumerate(hingeConfigurations)
    ΔM2min = [ΔM2peaks[c,i][1] for i in 1:length(ωRange)]
    ΔM2max = [ΔM2peaks[c,i][2] for i in 1:length(ωRange)]
    plot!(ωRange, ΔM2min, c=colors[c], lw=lw, label=false)
    plot!(ωRange, ΔM2max, c=colors[c], lw=lw, label=false)
    if config == :locked
        scatter!(ΔWRBM_locked_ref[1,:], ΔWRBM_locked_ref[2,:], c=colors[c], ms=ms, msw=msw, label=false)
    elseif config == :free
        scatter!(ΔWRBM_free_ref[1,:], ΔWRBM_free_ref[2,:], c=colors[c], ms=ms, msw=msw, label=false)
    end
end
display(plt_ΔWRBM)
savefig(string(absPath,"/HealyBaselineFFWTlockedVfreeOMCGustFloating_DWRBM.pdf"))

println("Finished HealyBaselineFFWTlockedVfreeOMCGustFloating.jl")
