using AeroBeams

# Beam
R,θ = 100, π/4
L = R*θ
A,Iy,Iz,J = 1, 1/12, 1/12, 1/6
E,G = 1e7, 5e6
nElem = 40
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E*A,GAy=G*A,GAz=G*A,GJ=G*J,EIy=E*Iy,EIz=E*Iz)],k=[0;0;1/R])

# BCs
F = 3000
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,conditions=Dict(:Ff3A=>F))

# Model
curvedCantileverStaticFollower = create_Model(name="curvedCantileverStaticFollower",beams=[beam],BCs=[clamp,tipForce],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options
σ0 = 0.0
σstep = 1e-2
NR = create_NewtonRaphson(initLoadFac=σ0,maxLoadFacStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=curvedCantileverStaticFollower,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
hist = problem.loadstepHistory
σVec = hist.σ
states = hist.states
tip_u1 = [states[i][nElem].nodal.u2[1] for i in 1:length(σVec)]
tip_u2 = [states[i][nElem].nodal.u2[2] for i in 1:length(σVec)]
tip_u3 = [states[i][nElem].nodal.u2[3] for i in 1:length(σVec)]

# Reference solution by Simo and Vu-Quoc (1986)
u1_ref = readdlm(pkgdir(AeroBeams,"test","referenceData","curvedCantileverStaticFollower", "u1.txt"))
u2_ref = readdlm(pkgdir(AeroBeams,"test","referenceData","curvedCantileverStaticFollower", "u2.txt"))
u3_ref = readdlm(pkgdir(AeroBeams,"test","referenceData","curvedCantileverStaticFollower", "u3.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/curvedCantileverStaticFollower"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,σ2plot=vcat(0.2:0.2:1),legendPos=(0.5,0.6),view=(135,15),save=true,savePath=string(relPath,"/curvedCantileverStaticFollower_deformation.pdf"))
display(deformationPlot)

# Plot configurations
colors = [:blue,:green,:orange]
labels = ["\$-u_1\$" "\$u_2\$" "\$u_3\$"]
ts = 10
fs = 16
lw = 2
ms = 3
msw = 0
gr()

# Plot normalized tip displacements over load steps
plt1 = plot(xlabel="\$F\$ [lb]", ylabel="Tip displacements [in]", xlims=[0,3000], xticks=collect(0:500:F), yticks=collect(-60:20:80), legend_position=(0.45, 0.4), tickfont=font(ts), guidefont=font(fs), legendfontsize=10)
plot!([NaN], [NaN], lc=:black,  lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Simo & Vu-Quoc (1986)")
for i=1:3
    plot!([NaN], [NaN], lc=colors[i], m=colors[i],  lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (y, c) in zip([-tip_u1, tip_u2, tip_u3], colors)
    plot!(σVec*F, y, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_ref[1,:],u2_ref[1,:],u3_ref[1,:]], [u1_ref[2,:],u2_ref[2,:],u3_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/curvedCantileverStaticFollower_disp.pdf"))

println("Finished curvedCantileverStaticFollower.jl")