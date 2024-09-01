using AeroBeams, ForwardDiff, DelimitedFiles

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)

# Mean airspeed and reduced frequency of airspeed oscillation
kᵤ = 0.2
Ma = 0.3
U₀ = Ma*atmosphere.a

# Wing surface
airfoil = create_Airfoil(name="flatPlate",Ma=Ma)
chord = 0.1
normSparPos = 1/4
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos,updateAirfoilParameters=false)

# Wing pitch
θ = 1*π/180

# Wing beam
L = 1
nElem = 1
∞ = 1e10
wing = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],rotationParametrization="E321",p0=[0;0;θ],aeroSurface=surf)

# BCs
clamp1 = create_BC(name="clamp1",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
timeVaryingFreestream = create_Model(name="timeVaryingFreestream",beams=[wing],BCs=[clamp1,clamp2],atmosphere=atmosphere)

# Set normalized airspeed amplitude range and initialize outputs
λᵤRange = collect(0.2:0.2:0.8)
t = Array{Vector{Float64}}(undef,length(λᵤRange))
tNorm = Array{Vector{Float64}}(undef,length(λᵤRange))
cn = Array{Vector{Float64}}(undef,length(λᵤRange))
cm = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot2 = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot3 = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot2Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot3Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
rangeLastCycle = Array{StepRange{Int64,Int64}}(undef,length(λᵤRange))

# Loop normalized airspeed amplitude
for (i,λᵤ) in enumerate(λᵤRange)
    # Update airspeed as a function of time
    ΔU = λᵤ*U₀
    ωᵤ = kᵤ*U₀/(chord/2)
    U = t -> U₀ + ΔU*sin.(ωᵤ*t)
    # Update velocity of basis A (and update model)
    set_motion_basis_A!(model=timeVaryingFreestream,v_A=t->[0;U(t);0])
    # Time variables
    T = 2π/ωᵤ
    cycles = 10
    tf = cycles*T
    Δt = T/200
    # Initial velocities update options
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,displayProgress=false, relaxFactor=0.5, Δt=Δt/1e3)
    # Create and solve problem
    global problem = create_DynamicProblem(model=timeVaryingFreestream,finalTime=tf,Δt=Δt,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
    solve!(problem)
    # Unpack numerical solution
    t[i] = problem.timeVector
    tNorm[i] = t[i]/T
    cn[i] = [problem.aeroVariablesOverTime[j][1].aeroCoefficients.cn for j in 1:length(t[i])]
    cm[i] = [problem.aeroVariablesOverTime[j][1].aeroCoefficients.cm for j in 1:length(t[i])]
    Vdot2[i] = [problem.elementalStatesRatesOverTime[j][1].Vdot[2] for j in 1:length(t[i])]
    Vdot3[i] = [problem.elementalStatesRatesOverTime[j][1].Vdot[3] for j in 1:length(t[i])]
    rangeLastCycle[i] = ceil(Int,(tf-T)/Δt):length(t[i])
    # Analytical solution for relative wind acceleration 
    Udot = t -> ForwardDiff.derivative(U,t)
    Vdot2Analytical[i] = Udot.(t[i]).*cos(θ)
    Vdot3Analytical[i] = -Udot.(t[i]).*sin(θ)
end

# Analytical solution for angular velocity and acceleration
HT = t -> AeroBeams.tangent_operator_transpose_WM([p(t);0;0])
HT_p = t -> AeroBeams.tangent_tensor_transpose_derivatives_extended_parameters([p(t);0;0])
HT_p1,HT_p2,HT_p3 = t -> HT_p(t)[1], t -> HT_p(t)[2], t -> HT_p(t)[3]
ΩAnalytical = t -> HT(t)*[pdot(t);0;0]
ΩdotAnalytical = t -> AeroBeams.mul3(HT_p1(t),HT_p2(t),HT_p3(t),[pdot(t);0;0])*[pdot(t);0;0] + HT(t)*[pddot(t);0;0]
Ω1Analytical = t -> ΩAnalytical(t)[1]
Ωdot1Analytical = t -> ΩdotAnalytical(t)[1]

# Quasi-steady normal force coefficient
cn_qs = t -> 2π*θ₀

# Load reference data
cnCFDLambda0_2 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "timeVaryingFreestream", "cnCFDLambda0_2.txt"))
cmCFDLambda0_2 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "timeVaryingFreestream", "cmCFDLambda0_2.txt"))
cnCFDLambda0_4 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "timeVaryingFreestream", "cnCFDLambda0_4.txt"))
cmCFDLambda0_4 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "timeVaryingFreestream", "cmCFDLambda0_4.txt"))
cnCFDLambda0_6 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "timeVaryingFreestream", "cnCFDLambda0_6.txt"))
cmCFDLambda0_6 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "timeVaryingFreestream", "cmCFDLambda0_2.txt"))
cnCFDLambda0_8 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "timeVaryingFreestream", "cnCFDLambda0_8.txt"))
cmCFDLambda0_8 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "timeVaryingFreestream", "cmCFDLambda0_8.txt"))

cnCFD = Array{Matrix{Float64}}(undef,4)
cmCFD = Array{Matrix{Float64}}(undef,4)
cnCFD[1] = cnCFDLambda0_2
cmCFD[1] = cmCFDLambda0_2
cnCFD[2] = cnCFDLambda0_4
cmCFD[2] = cmCFDLambda0_4
cnCFD[3] = cnCFDLambda0_6
cmCFD[3] = cmCFDLambda0_6
cnCFD[4] = cnCFDLambda0_8
cmCFD[4] = cmCFDLambda0_8

println("Finished timeVaryingFreestream.jl")