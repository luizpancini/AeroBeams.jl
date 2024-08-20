using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# This problem was defined by Argyris et al. - Finite Element Method: the natural approach (1979)

# Beam frame
L,b,H = 2.4e-1,6e-4,3e-2
A,Iy,Iz = b*H,b*H^3/12,H*b^3/12
J = Iy+Iz
Ksz,Ksy,Kt = 5/6,1/20,1/50
E = 7.124e10
ν = 0.31
G = E/(2*(1+ν))
EA,GAy,GAz,GJ,EIy,EIz = E*A,G*A*Ksy,G*A*Ksz,G*J*Kt,E*Iy,E*Iz
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
nElem = 10
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,C=[stiffnessMatrix])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,C=[stiffnessMatrix],rotationParametrization="E321",p0=[0;π/2;0])

# BCs
δ = 1e-3
F = 2
clamp = create_BC(name="clamp",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipMisalignedForce = create_BC(name="tipMisalignedForce",beam=beam2,node=nElem+1,types=["F1A","F2A"],values=[F,F*δ])

# Model
rightAngledFrameBuckling = create_Model(name="rightAngledFrameBuckling",beams=[beam1,beam2],BCs=[clamp,tipMisalignedForce])

# Set system solver options
σ0 = 0.0
σstep = 0.01
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=rightAngledFrameBuckling,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u2 = [problem.nodalStatesOverσ[i][end].u_n2[2] for i in 1:length(σVector)]

# Plot deformed shape
relPath = "/test/outputs/figures/rightAngledFrameBuckling"
absPath = string(pwd(),relPath)
mkpath(absPath)
deformationPlot = plot_steady_deformation(problem,save=true,view=(45,20),savePath=string(relPath,"/rightAngledFrameBuckling_deformation.pdf"))
display(deformationPlot)

# Plot normalized tip out-of-plane displacement over load steps
gr()
plt1 = plot(ylabel="\$F\$ [N]", xlabel="Tip \$u_2/L\$", title="Tip out-of-plane displacement")
plot!(tip_u2/L, σVector*F, lw=2, label=false)
display(plt1)
savefig(string(absPath,"/rightAngledFrameBuckling_disp.pdf"))

println("Finished rightAngledFrameBuckling.jl")