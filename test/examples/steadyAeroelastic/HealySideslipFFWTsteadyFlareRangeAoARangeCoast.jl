using AeroBeams

# Flare angle range [rad]
ΛRange = deg2rad.([10,20,30])

# Root pitch angle range [rad]
θRange = deg2rad.(vcat(-20:1:30))

# Airspeed [m/s]
U = 22

# Discretization
nElementsInner = 15
nElementsFFWT = 5

# Tip loss options
hasTipCorrection = true
tipLossDecayFactor = 10

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Initialize model
model = create_HealySideslipFFWT(flareAngle=0,pitchAngle=0,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)

# Fixed properties of the model
nElem = nElementsInner + nElementsFFWT
beam = model.beams[1]
r_n1 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[1], 
                    beam.elements[e].assemblyData.r_n2[1]] for e in 1:nElem)
r_n2 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[2], 
                    beam.elements[e].assemblyData.r_n2[2]] for e in 1:nElem)
r_n3 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[3], 
                    beam.elements[e].assemblyData.r_n2[3]] for e in 1:nElem)
x1_n = reduce(vcat,[beam.elements[e].geometry.x1_n1[1],
                    beam.elements[e].geometry.x1_n2[1]] for e in 1:nElem)
x1_e = [beam.elements[e].geometry.x1 for e in 1:nElem]

# Initialize outputs
u1 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
u2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
u3 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
p1 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
p2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
p3 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
M2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
α = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
cn = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
pHinge = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange))
ϕHinge = Array{Float64}(undef,length(ΛRange),length(θRange))

# Loop flare angle
for (i,Λ) in enumerate(ΛRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Loop root pitch angle
    for (j,θ) in enumerate(θRange)
        println("Solving for Λ = $(round(Int,rad2deg(Λ))) deg, θ = $(round(Int,rad2deg(θ))) deg")
        # Update model
        model = create_HealySideslipFFWT(flareAngle=Λ,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)
        # Create and solve problem
        problem = create_SteadyProblem(model=model,systemSolver=systemSolver,x0=x0)
        solve!(problem)
        # Get outputs
        state = problem.loadstepHistory.states[end]
        hacState = problem.loadstepHistory.hacStates[end]
        aeroVars = problem.loadstepHistory.aeroVariables[end]
        u1[i,j] = reduce(vcat,[state[e].nodal.u1[1],
                            state[e].nodal.u2[1]] for e in 1:nElem)
        u2[i,j] = reduce(vcat,[state[e].nodal.u1[2],
                            state[e].nodal.u2[2]] for e in 1:nElem)
        u3[i,j] = reduce(vcat,[state[e].nodal.u1[3],
                            state[e].nodal.u2[3]] for e in 1:nElem)
        p1[i,j] = reduce(vcat,[state[e].nodal.p1[1],
                            state[e].nodal.p2[1]] for e in 1:nElem)
        p2[i,j] = reduce(vcat,[state[e].nodal.p1[2],
                            state[e].nodal.p2[2]] for e in 1:nElem)
        p3[i,j] = reduce(vcat,[state[e].nodal.p1[3],
                            state[e].nodal.p2[3]] for e in 1:nElem)
        M2[i,j] = reduce(vcat,[state[e].nodal.M1[2],
                            state[e].nodal.M2[2]] for e in 1:nElem)
        α[i,j] = [rad2deg(aeroVars[e].unsteadyKin.αₑ) for e in 1:nElem]
        cn[i,j] = [aeroVars[e].aeroCoefficients.cn for e in 1:nElem]
        pHinge[i,j] = hacState[1].pH
        ϕHinge[i,j] = rad2deg(hacState[1].ϕ)
        # Set initial guess for next iteration
        x0 = problem.systemArrays.x
    end
end

# Load reference data
flare10_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare10_exp.txt")
flare10_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare10_num.txt")
flare20_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare20_exp.txt")
flare20_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare20_num.txt")
flare30_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare30_exp.txt")
flare30_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/flare30_num.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/HealySideslipFFWTsteadyFlareRangeAoARangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(ΛRange), categorical=true)
ts = 12
fs = 16
lfs = 12
lw = 2
ms = 6
msw = 0

# Coast angle vs pitch angle for each flare angle
plt = plot(xlabel="Root pitch angle [deg]", ylabel="Coast angle [deg]", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, xlims=[-20,30], ylims=[-90,90], yticks=-90:30:90, legend=:topleft)
scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Healy (2023) - Exp.")
plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (i,Λ) in enumerate(ΛRange)
    plot!(rad2deg.(θRange), -ϕHinge[i,:], lw=lw, ls=:solid, c=colors[i], label="\$\\Lambda = $(round(Int,rad2deg(Λ))) \\degree\$")
    if i==1
        plot!(flare10_num[1,:], flare10_num[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
        scatter!(flare10_exp[1,:], flare10_exp[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif i==2
        plot!(flare20_num[1,:], flare20_num[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
        scatter!(flare20_exp[1,:], flare20_exp[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    else
        plot!(flare30_num[1,:], flare30_num[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
        scatter!(flare30_exp[1,:], flare30_exp[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt)
savefig(string(absPath,"/HealySideslipFFWTsteadyFlareRangeAoARangeCoast_coastAngle.pdf"))

println("Finished HealySideslipFFWTsteadyFlareRangeAoARangeCoast.jl")