using AeroBeams, LinearAlgebra, Plots, ColorSchemes, BenchmarkTools

## User inputs (problem definition)
#-------------------------------------------------------------------------------
# Geometric properties
L = 1
b1,b2,h = 5e-2,15e-2,5e-2
A1,Iy1,Iz1 = b1*h,b1*h^3/12,h*b1^3/12
A2,Iy2,Iz2 = b2*h,b2*h^3/12,h*b2^3/12
J1 = Is1 = A1^4/(Iy1+Iz1)/40
J2 = Is2 = A2^4/(Iy2+Iz2)/40
Kt1,Kt2 = 10/9,3.0864
α = [0; 0; π/2; 0; 0; π/2; 0; 0; 0; 0; 0; 0; π/2; 0; 0; π/2]
β = -[π/2; π/2; 0; -π/2; -π/2; 0; 0; 0; 0; 0; -π/2; π/2; 0; -π/2; -π/2; 0]

# Material properties
E = 219.9e9
ν = 0.25
G = E/(2*(1+ν))
ρ = 7.9e3
∞ = 1e12

# Stiffness and inertia matrices
stiffnessMatrices = [diagm([E*A1,∞,∞,G*J1*Kt1,E*Iy1,E*Iz1]),diagm([E*A2,∞,∞,G*J2*Kt2,E*Iy2,E*Iz2])]
inertiaMatrices = [diagm([ρ*A1,ρ*A1,ρ*A1,ρ*Is1,ρ*Iy1,ρ*Iz1]),diagm([ρ*A2,ρ*A2,ρ*A2,ρ*Is2,ρ*Iy2,ρ*Iz2])]

# Number of elements for each beam
nElem = 10

## Problem setup
#-------------------------------------------------------------------------------
# Beams
beams = Vector{Beam}(undef,16)

beams[1] = create_Beam(name="beam1",length=L,nElements=nElem,C=[stiffnessMatrices[1]],I=[inertiaMatrices[1]],rotationParametrization="E321",p0=[α[1];β[1];0.0]) 

beams[2] = create_Beam(name="beam2",length=L,nElements=nElem,C=[stiffnessMatrices[1]],I=[inertiaMatrices[1]],rotationParametrization="E321",p0=[α[2];β[2];0.0],connectedBeams=[beams[1]],connectedNodesThis=[1],connectedNodesOther=[nElem+1])

beams[3] = create_Beam(name="beam3",length=L,nElements=nElem,C=[stiffnessMatrices[2]],I=[inertiaMatrices[2]],rotationParametrization="E321",p0=[α[3];β[3];0.0],connectedBeams=[beams[2]],connectedNodesThis=[1],connectedNodesOther=[nElem+1])

beams[4] = create_Beam(name="beam4",length=L,nElements=nElem,C=[stiffnessMatrices[1]],I=[inertiaMatrices[1]],rotationParametrization="E321",p0=[α[4];β[4];0.0],connectedBeams=[beams[3]],connectedNodesThis=[1],connectedNodesOther=[nElem+1])

beams[5] = create_Beam(name="beam5",length=L,nElements=nElem,C=[stiffnessMatrices[1]],I=[inertiaMatrices[1]],rotationParametrization="E321",p0=[α[5];β[5];0.0],connectedBeams=[beams[4]],connectedNodesThis=[1],connectedNodesOther=[nElem+1])

beams[6] = create_Beam(name="beam6",length=L,nElements=nElem,C=[stiffnessMatrices[2]],I=[inertiaMatrices[2]],rotationParametrization="E321",p0=[α[6];β[6];0.0],connectedBeams=[beams[1],beams[4]],connectedNodesThis=[1,nElem+1],connectedNodesOther=[nElem+1,nElem+1])

beams[7] = create_Beam(name="beam7",length=L,nElements=nElem,C=[stiffnessMatrices[2]],I=[inertiaMatrices[2]],rotationParametrization="E321",p0=[α[7];β[7];0.0],connectedBeams=[beams[1]],connectedNodesThis=[1],connectedNodesOther=[nElem+1])

beams[8] = create_Beam(name="beam8",length=L,nElements=nElem,C=[stiffnessMatrices[2]],I=[inertiaMatrices[2]],rotationParametrization="E321",p0=[α[8];β[8];0.0],connectedBeams=[beams[4]],connectedNodesThis=[1],connectedNodesOther=[nElem+1])

beams[9] = create_Beam(name="beam9",length=L,nElements=nElem,C=[stiffnessMatrices[2]],I=[inertiaMatrices[2]],rotationParametrization="E321",p0=[α[9];β[9];0.0],connectedBeams=[beams[2]],connectedNodesThis=[1],connectedNodesOther=[nElem+1])

beams[10] = create_Beam(name="beam10",length=L,nElements=nElem,C=[stiffnessMatrices[2]],I=[inertiaMatrices[2]],rotationParametrization="E321",p0=[α[10];β[10];0.0],connectedBeams=[beams[3]],connectedNodesThis=[1],connectedNodesOther=[nElem+1])

beams[11] = create_Beam(name="beam11",length=L,nElements=nElem,C=[stiffnessMatrices[1]],I=[inertiaMatrices[1]],rotationParametrization="E321",p0=[α[11];β[11];0.0],connectedBeams=[beams[7]],connectedNodesThis=[1],connectedNodesOther=[nElem+1])

beams[12] = create_Beam(name="beam12",length=L,nElements=nElem,C=[stiffnessMatrices[1]],I=[inertiaMatrices[1]],rotationParametrization="E321",p0=[α[12];β[12];0.0],connectedBeams=[beams[7],beams[9]],connectedNodesThis=[1,nElem+1],connectedNodesOther=[nElem+1,nElem+1])

beams[13] = create_Beam(name="beam13",length=L,nElements=nElem,C=[stiffnessMatrices[2]],I=[inertiaMatrices[2]],rotationParametrization="E321",p0=[α[13];β[13];0.0],connectedBeams=[beams[9],beams[10]],connectedNodesThis=[1,nElem+1],connectedNodesOther=[nElem+1,nElem+1])

beams[14] = create_Beam(name="beam14",length=L,nElements=nElem,C=[stiffnessMatrices[1]],I=[inertiaMatrices[1]],rotationParametrization="E321",p0=[α[14];β[14];0.0],connectedBeams=[beams[13],beams[8]],connectedNodesThis=[1,nElem+1],connectedNodesOther=[nElem+1,nElem+1])

beams[15] = create_Beam(name="beam15",length=L,nElements=nElem,C=[stiffnessMatrices[1]],I=[inertiaMatrices[1]],rotationParametrization="E321",p0=[α[15];β[15];0.0],connectedBeams=[beams[8]],connectedNodesThis=[1],connectedNodesOther=[nElem+1])

beams[16] = create_Beam(name="beam16",length=L,nElements=nElem,C=[stiffnessMatrices[2]],I=[inertiaMatrices[2]],rotationParametrization="E321",p0=[α[16];β[16];0.0],connectedBeams=[beams[7],beams[8]],connectedNodesThis=[1,nElem+1],connectedNodesOther=[nElem+1,nElem+1])

# BCs
clamps = Vector{BC}(undef,4)
clamps[1] = create_BC(name="clamp1",beam=beams[1],node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamps[2] = create_BC(name="clamp2",beam=beams[5],node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamps[3] = create_BC(name="clamp3",beam=beams[11],node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamps[4] = create_BC(name="clamp4",beam=beams[15],node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Create model
twoStoryFrame = create_Model(name="twoStoryFrame",beams=beams,BCs=clamps,units=create_UnitsSystem(length="m",frequency="Hz"))

# Create and solve eigenproblem
problem = create_EigenProblem(model=twoStoryFrame,nModes=4,getLinearSolution=true)
solve!(problem)

# Get frequencies
freqs = problem.frequenciesOscillatory

# Reference frequencies (in Hz) by PETYT - Introduction to Finite Element Vibration Analysis - [2nd Ed.] (2010)
refFreqs = [11.8; 34.1]

# Display relative errors
ϵ_rel = freqs[[1,4]]/(2π)./refFreqs .- 1.0
println("Relative errors: $ϵ_rel")

# Plot mode shapes
relPath = "/test/outputs/figures/twoStoryFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)

modesPlot = plot_mode_shapes(problem,scale=1,view=(45,30),legendPos=(0.3,0.1),frequencyLabel="frequency",save=true,savePath=string(relPath,"/twoStoryFrame_modeShapes.pdf"))
display(modesPlot)

println("Finished twoStoryFrame.jl")