# Script to replicate the problem defined by CHENG et al. (2025) - Nonlinear Multibody Modelling of Flexible Aircraft with Flared Hinged Wing Tips [JoA]

using AeroBeams, DelimitedFiles

# Slant (flare) angle configurations
ΛRange = π/2*[0; 1/2; 1]

# Solution method for hinge axis constraint
solutionMethod = "addedResidual"
updateAllDOFinResidual = true

# Gravity
g = 9.80665

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
tf = 2
Δt = tf/5e2

# System solver
σ0 = 1
maxIter = 100
absTol = 1e-8
relTol = 1e-7
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,absoluteTolerance=absTol,relativeTolerance=relTol)

# Initialize outputs
problem = Array{DynamicProblem}(undef,length(ΛRange))
t = Array{Vector{Float64}}(undef,length(ΛRange))
u1_tip = Array{Vector{Float64}}(undef,length(ΛRange))
u2_tip = Array{Vector{Float64}}(undef,length(ΛRange))
u3_tip = Array{Vector{Float64}}(undef,length(ΛRange))

# Sweep flare angle
for (i,Λ) in enumerate(ΛRange)
    println("Solving for Λ = $(round(Int,Λ*180/π)) deg")
    # Hinge node DOFs
    if Λ == 0
        hingedNodesDoF = [false,true,false]
    elseif Λ == π/2
        hingedNodesDoF = falses(3)
    else
        hingedNodesDoF = trues(3)
    end
    # Beam
    beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA,GAy=GA,GAz=GA,EIy=EI,EIz=EI,GJ=GJ)],I=[inertia_matrix(ρA=ρ*A,ρIy=ρ*Iy,ρIz=ρ*Iz)],hingedNodes=[hingeNode],hingedNodesDoF=[hingedNodesDoF])
    # Hinge axis constraint
    hingeAxisConstraints = Vector{HingeAxisConstraint}()
    if 0 < Λ < π/2
        localHingeAxis = rotation_tensor_E321([-Λ; 0; 0]) * AeroBeams.a2
        push!(hingeAxisConstraints,create_HingeAxisConstraint(solutionMethod=solutionMethod,updateAllDOFinResidual=updateAllDOFinResidual,beam=beam,localHingeAxis=localHingeAxis))
    end
    # Spring around hinge
    kT = 1e-4
    kOOP = 1e-4
    kIP = 1e-4
    spring = create_Spring(elementsIDs=[inboardElem,outboardElem],nodesSides=[1,2],kp=[kT,kOOP,kIP])
    if 0 < Λ < π/2
        add_spring_to_beams!(beams=[beam,beam],spring=spring)
    end
    # BCs
    support = create_BC(name="support",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
    # Model
    slantedDoublePendulum = create_Model(name="slantedDoublePendulum",beams=[beam],BCs=[support],gravityVector=[0,0,-g],hingeAxisConstraints=hingeAxisConstraints)
    # Create and solve the problem
    problem[i] = create_DynamicProblem(model=slantedDoublePendulum,finalTime=tf,Δt=Δt,systemSolver=NR,skipInitialStatesUpdate=true)
    solve!(problem[i])
    # Unpack numerical solution
    t[i] = problem[i].savedTimeVector
    u1_tip[i] = [problem[i].nodalStatesOverTime[j][end].u_n2[1] for j in eachindex(t[i])]
    u2_tip[i] = [problem[i].nodalStatesOverTime[j][end].u_n2[2] for j in eachindex(t[i])]
    u3_tip[i] = [problem[i].nodalStatesOverTime[j][end].u_n2[3] for j in eachindex(t[i])]
end

# Load reference data
Λ0_r1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda0_r1.txt")
Λ0_r3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda0_r3.txt")

Λ45_r1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda45_r1.txt")
Λ45_r2_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda45_r2.txt")
Λ45_r3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda45_r3.txt")

Λ90_r1_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda90_r1.txt")
Λ90_r3_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/ChengSlantedDoublePendulum/Lambda90_r3.txt")

println("Finished slantedDoublePendulum.jl")