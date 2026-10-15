using AeroBeams

# Flare angle range [rad]
ΛRange = deg2rad.([10,20,30])

# Root pitch angle range
θRange = deg2rad.([0,5,10])

# Airspeed range [m/s]
URange = collect(10:1:35)

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
u1 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
u2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
u3 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
p1 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
p2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
p3 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
M2 = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
α = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
cn = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
pHinge = Array{Vector{Float64}}(undef,length(ΛRange),length(θRange),length(URange))
ϕHinge = Array{Float64}(undef,length(ΛRange),length(θRange),length(URange))

# Loop flare angle
for (i,Λ) in enumerate(ΛRange)
    # Loop root pitch angle
    for (j,θ) in enumerate(θRange)
        # Initialize guess solution
        x0 = zeros(0)
        # Loop airspeed
        for (k,U) in enumerate(URange)
            println("Solving for Λ = $(round(Int,rad2deg(Λ))) deg, θ = $(round(Int,rad2deg(θ))) deg, U = $U m/s")
            # Update model
            model = create_HealySideslipFFWT(flareAngle=Λ,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)
            # Create and solve problem
            problem = create_SteadyProblem(model=model,systemSolver=systemSolver,x0=x0)
            solve!(problem)
            # Get outputs
            state = problem.loadstepHistory.states[end]
            hacState = problem.loadstepHistory.hacStates[end]
            aeroVars = problem.loadstepHistory.aeroVariables[end]
            u1[i,j,k] = reduce(vcat,[state[e].nodal.u1[1],
                              state[e].nodal.u2[1]] for e in 1:nElem)
            u2[i,j,k] = reduce(vcat,[state[e].nodal.u1[2],
                              state[e].nodal.u2[2]] for e in 1:nElem)
            u3[i,j,k] = reduce(vcat,[state[e].nodal.u1[3],
                              state[e].nodal.u2[3]] for e in 1:nElem)
            p1[i,j,k] = reduce(vcat,[state[e].nodal.p1[1],
                              state[e].nodal.p2[1]] for e in 1:nElem)
            p2[i,j,k] = reduce(vcat,[state[e].nodal.p1[2],
                              state[e].nodal.p2[2]] for e in 1:nElem)
            p3[i,j,k] = reduce(vcat,[state[e].nodal.p1[3],
                              state[e].nodal.p2[3]] for e in 1:nElem)
            M2[i,j,k] = reduce(vcat,[state[e].nodal.M1[2],
                              state[e].nodal.M2[2]] for e in 1:nElem)
            α[i,j,k] = [rad2deg(aeroVars[e].unsteadyKin.αₑ) for e in 1:nElem]
            cn[i,j,k] = [aeroVars[e].aeroCoefficients.cn for e in 1:nElem]
            pHinge[i,j,k] = hacState[1].pH
            ϕHinge[i,j,k] = rad2deg(hacState[1].ϕ)
        end
    end
end

# Load reference data
flare10_aoa0_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa0_exp.txt")
flare10_aoa0_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa0_num.txt")
flare10_aoa5_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa5_exp.txt")
flare10_aoa5_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa5_num.txt")
flare10_aoa10_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa10_exp.txt")
flare10_aoa10_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare10_aoa10_num.txt")
flare20_aoa0_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa0_exp.txt")
flare20_aoa0_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa0_num.txt")
flare20_aoa5_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa5_exp.txt")
flare20_aoa5_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa5_num.txt")
flare20_aoa10_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa10_exp.txt")
flare20_aoa10_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare20_aoa10_num.txt")
flare30_aoa0_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa0_exp.txt")
flare30_aoa0_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa0_num.txt")
flare30_aoa5_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa5_exp.txt")
flare30_aoa5_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa5_num.txt")
flare30_aoa10_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa10_exp.txt")
flare30_aoa10_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/flare30_aoa10_num.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
ts = 12
fs = 16
lfs = 10
lw = 2
ms = 6
msw = 0

# Coast angle vs airspeed for each AoA, for several flare angles
for (i,Λ) in enumerate(ΛRange)
    if i==1
        figurePath = "/HealySideslipFFWTsteadyURangeAoARangeCoast_phi_flare10.pdf"
        aoa0_exp = flare10_aoa0_exp
        aoa0_num = flare10_aoa0_num
        aoa5_exp = flare10_aoa5_exp
        aoa5_num = flare10_aoa5_num
        aoa10_exp = flare10_aoa10_exp
        aoa10_num = flare10_aoa10_num
    elseif i==2
        figurePath = "/HealySideslipFFWTsteadyURangeAoARangeCoast_phi_flare20.pdf"
        aoa0_exp = flare20_aoa0_exp
        aoa0_num = flare20_aoa0_num
        aoa5_exp = flare20_aoa5_exp
        aoa5_num = flare20_aoa5_num
        aoa10_exp = flare20_aoa10_exp
        aoa10_num = flare20_aoa10_num
    else
        figurePath = "/HealySideslipFFWTsteadyURangeAoARangeCoast_phi_flare30.pdf"
        aoa0_exp = flare30_aoa0_exp
        aoa0_num = flare30_aoa0_num
        aoa5_exp = flare30_aoa5_exp
        aoa5_num = flare30_aoa5_num
        aoa10_exp = flare30_aoa10_exp
        aoa10_num = flare30_aoa10_num
    end
    plt = plot(xlabel="Airspeed [m/s]", ylabel="Coast angle [deg]", title="Flare angle: $(round(Int,rad2deg(Λ)))°", xlims=[10,35], ylims=[-90,45], yticks=-90:30:45, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:bottomright)
    scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Healy (2023) - Exp.")
    plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
    plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
    for (j,θ) in enumerate(θRange)
        plot!(URange, -ϕHinge[i,j,:], lw=lw, ls=:solid, c=colors[j], label="\$\\theta = $(round(Int,rad2deg(θ))) \\degree\$")
        if j==1
            plot!(aoa0_num[1,:], aoa0_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoa0_exp[1,:], aoa0_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)
        elseif j==2
            plot!(aoa5_num[1,:], aoa5_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoa5_exp[1,:], aoa5_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)
        else
            plot!(aoa10_num[1,:], aoa10_num[2,:], lw=lw, ls=:dash, c=colors[j], label=false)
            scatter!(aoa10_exp[1,:], aoa10_exp[2,:], ms=ms, msw=msw, c=colors[j], label=false)
        end
    end
    display(plt)
    savefig(string(absPath,figurePath))
end

println("Finished HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast.jl")