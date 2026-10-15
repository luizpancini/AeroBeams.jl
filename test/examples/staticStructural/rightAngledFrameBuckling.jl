using AeroBeams

# This problem was defined by Argyris et al. - Finite Element Method: the natural approach (1979)

# Beam frame
L,b,H = 2.4e-1, 6e-4, 3e-2
A,Iy,Iz = b*H, b*H^3/12, H*b^3/12
J = Iy+Iz
Ksz,Ksy,Kt = 5/6, 1/20, 1/50
E = 7.124e10
ν = 0.31
G = E/(2*(1+ν))
S = isotropic_stiffness_matrix(EA=E*A,GAy=G*A*Ksy,GAz=G*A*Ksz,GJ=G*J*Kt,EIy=E*Iy,EIz=E*Iz)
nElem = 10
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,S=[S])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,S=[S],p0=transform_to_WM([0;π/2;0],from=:E321))

# BCs
δ = 1e-3
F = 2
clamp = create_BC(name="clamp",beam=beam1,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipMisalignedForce = create_BC(name="tipMisalignedForce",beam=beam2,node=nElem+1,conditions=Dict(:F1A=>F,:F2A=>F*δ))

# Model
rightAngledFrameBuckling = create_Model(name="rightAngledFrameBuckling",beams=[beam1,beam2],BCs=[clamp,tipMisalignedForce])

# Set system solver options
σ0 = 0.0
σstep = 0.01
NR = create_NewtonRaphson(initLoadFac=σ0,maxLoadFacStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=rightAngledFrameBuckling,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
hist = problem.loadstepHistory
σVec = hist.σ
states = hist.states
tip_u2 = [states[i][end].nodal.u2[2] for i in 1:length(σVec)]

## Plots
# Set paths
relPath = "/test/outputs/figures/rightAngledFrameBuckling"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,view=(45,20),savePath=string(relPath,"/rightAngledFrameBuckling_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Plot normalized tip out-of-plane displacement over load steps
plt1 = plot(ylabel="\$F\$ [kN]", xlabel="Tip \$u_2/L\$", title="Tip out-of-plane displacement")
plot!(tip_u2/L, σVec*F/1e3, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/rightAngledFrameBuckling_disp.pdf"))

println("Finished rightAngledFrameBuckling.jl")