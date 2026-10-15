# Script to replicate the problem defined by CHENG et al. (2025) - Nonlinear Multibody Modelling of Flexible Aircraft with Flared Hinged Wing Tips [JoA]
using AeroBeams

# Slant (flare) angle configurations [rad]
ΛRange = π/2*[0, 1/2, 1]

# Discretization
nElem = 10

# Dimensions
L = 1
d = 2e-2
A = d^2
Iy = Iz = d^4/12

# Sectional stiffness and inertia properties
EA = 2.8e7
GA = 1.037e7
EI = 9.333e2
GJ = 6.914e2
ρ = 2.7e3

# Slanted (flared) hinge
hingeNode = div(nElem,2)+1
inboardElem = hingeNode - 1
outboardElem = inboardElem + 1

# Time variables
nSteps = 1000
tf = 2
Δt = tf/nSteps

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true,displayStatus=false)

# Initialize outputs
problem = Array{DynamicProblem}(undef,length(ΛRange))
t = Array{Vector{Float64}}(undef,length(ΛRange))
u1_tip = Array{Vector{Float64}}(undef,length(ΛRange))
u2_tip = Array{Vector{Float64}}(undef,length(ΛRange))
u3_tip = Array{Vector{Float64}}(undef,length(ΛRange))

# Sweep flare angle
for (i,Λ) in enumerate(ΛRange)
    println("Solving for Λ = $(round(Int,rad2deg(Λ))) deg")
    # Hinge node DOFs
    hingedNodesDoF = if Λ == 0
        [false,true,false]
    elseif Λ == π/2
        falses(3)
    else
        trues(3)
    end
    # Beam
    beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA,GAy=GA,GAz=GA,EIy=EI,EIz=EI,GJ=GJ)],I=[inertia_matrix(ρA=ρ*A,ρIy=ρ*Iy,ρIz=ρ*Iz)],hingedNodes=[hingeNode],hingedNodesDoF=[hingedNodesDoF])
    # Hinge axis constraint
    hingeAxisConstraints = Vector{HingeAxisConstraint}()
    if 0 < Λ < π/2
        localHingeAxis = rotation_tensor_E321([-Λ, 0, 0]) * AeroBeams.a2
        push!(hingeAxisConstraints,create_HingeAxisConstraint(beam=beam,localHingeAxis=localHingeAxis))
    end
    # BCs
    support = create_BC(name="support",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
    # Model
    slantedDoublePendulum = create_Model(name="slantedDoublePendulum",beams=[beam],BCs=[support],gravityVector=[0,0,-standard_atmosphere(0).g],hingeAxisConstraints=hingeAxisConstraints)
    # Create and solve the problem
    problem[i] = create_DynamicProblem(model=slantedDoublePendulum,finalTime=tf,Δt=Δt,systemSolver=systemSolver,skipInitialStatesUpdate=true)
    solve!(problem[i])
    # @profview_allocs solve!(problem[i]) sample_rate=2e-2
    # @time solve!(problem[i])
    # Unpack numerical solution
    sol = problem[i].solutionHistory
    states = sol.states
    t[i] = sol.t
    u1_tip[i] = [states[j][end].nodal.u2[1] for j in eachindex(t[i])]
    u2_tip[i] = [states[j][end].nodal.u2[2] for j in eachindex(t[i])]
    u3_tip[i] = [states[j][end].nodal.u2[3] for j in eachindex(t[i])]
end

# Load reference data
Λ0_r1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda0_r1.txt")
Λ0_r3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda0_r3.txt")

Λ45_r1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda45_r1.txt")
Λ45_r2_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda45_r2.txt")
Λ45_r3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda45_r3.txt")

Λ90_r1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda90_r1.txt")
Λ90_r3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda90_r3.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/ChengSlantedDoublePendulum"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 4
msw = 0
ts = 12
fs = 16
lfs = 10
Λlabels = ["\$\\Lambda = $(round(Int,rad2deg(Λ)))^\\circ\$" for Λ in ΛRange]
colors = cgrad(:rainbow, length(ΛRange), categorical=true)
gr()

# Normalized tip u1 displacement
plt_u1 = plot(xlabel="\$t\$ [s]", ylabel="\$u_1/L\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:top)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, msw=msw, label="Cheng et al. (2025)")
for (i,Λ) in enumerate(ΛRange)
    plot!(t[i], u1_tip[i]/L, lw=lw, c=colors[i], label=Λlabels[i])
    if Λ == 0
        scatter!(Λ0_r1_ref[1,:], Λ0_r1_ref[2,:].-L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif Λ == π/4
        scatter!(Λ45_r1_ref[1,:], Λ45_r1_ref[2,:].-L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif Λ == π/2
        scatter!(Λ90_r1_ref[1,:], Λ90_r1_ref[2,:].-L, ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_u1)
savefig(plt_u1, string(absPath,"/ChengSlantedDoublePendulum_u1.pdf"))

# Normalized tip u2 displacement
plt_u2 = plot(xlabel="\$t\$ [s]", ylabel="\$u_2/L\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (i,Λ) in enumerate(ΛRange)
    plot!(t[i], u2_tip[i]/L, lw=lw, c=colors[i], label=false)
    if Λ == π/4
        scatter!(Λ45_r2_ref[1,:], Λ45_r2_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_u2)
savefig(plt_u2, string(absPath,"/ChengSlantedDoublePendulum_u2.pdf"))

# Normalized tip u3 displacement
plt_u3 = plot(palette=colors, xlabel="\$t\$ [s]", ylabel="\$u_3/L\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (i,Λ) in enumerate(ΛRange)
    plot!(t[i], u3_tip[i]/L, lw=lw, c=colors[i], label=false)
    if Λ == 0
        scatter!(Λ0_r3_ref[1,:], Λ0_r3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif Λ == π/4
        scatter!(Λ45_r3_ref[1,:], Λ45_r3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif Λ == π/2
        scatter!(Λ90_r3_ref[1,:], Λ90_r3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_u3)
savefig(plt_u3, string(absPath,"/ChengSlantedDoublePendulum_u3.pdf"))

println("Finished ChengSlantedDoublePendulum.jl")