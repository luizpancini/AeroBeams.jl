using AeroBeams, ForwardDiff

# Reference data source ("NASA" or "GU")
source = "NASA"

# Frame
frame = 10_022

# Number of cycles
nCycles = 4

# Starting position of the cycle (:begin or :middle)
startingPosition = :middle

# Flag to apply a transition function to pitch profile
applyTransitionFunction = true

# Aerodynamic solver
aeroSolver = BLi(circulatoryIndicialFunction=:Beddoes,incompressibleInertialLoads=false)

# Frame data
if source == "NASA"
    airfoil,a₀,a₁,b,k,Ma,U = NASA_frames_loader(frame)
elseif source == "GU"
    airfoil,a₀,a₁,b,k,Ma,U = GU_frames_loader(frame)
else
    error("Unavailable source")
end

# Dimensional frequency and cycle time
ω = k*U/b
τ = 2π/ω

# Pitch profile
θ = startingPosition == :begin ? 
    t -> a₁*(1+sin(ω*(t-τ/4))) : 
    t -> a₁*sin(ω*t)
if applyTransitionFunction
    θ = let θ_old = θ
        t -> θ_old(t) * (t/τ)^10 / (1 + (t/τ)^10)
    end
end
p = t -> 4*tan(θ(t)/4)
pdot = t -> ForwardDiff.derivative(p,t)

# Aerodynamic surface
surf = create_AeroSurface(ppSolver=aeroSolver,airfoil=airfoil,c=2*b,normSparPos=0.25,updateAirfoilParameters=false)

# Wing (rigid)
L = 1
nElem = 1
∞ = 1e10
ρA,ρI = 10, 1e-2
p0₃ = startingPosition == :begin ? a₀-a₁ : a₀
wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=ρA,ρIy=ρI,ρIz=ρI)],p0=transform_to_WM([0,0,p0₃],from=:E321),surface=surf)

# Set initial conditions
add_initial_conditions_to_beam!(beam=wing,conditions=Dict(:pdot0_of_x1=>x1->[pdot(0), 0, 0]))

# BCs
driver = create_BC(name="driver",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>t -> p(t),:p2A=>0,:p3A=>0))
journal = create_BC(name="journal",beam=wing,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p2A=>0,:p3A=>0))

# Model
pitchingAirfoilDSModelTest = create_Model(name="pitchingAirfoilDSModelTest",beams=[wing],BCs=[driver,journal],v_A=[0,U,0])

# Time variables
nTimeStepsCycle = 500
Δt = τ / nTimeStepsCycle
tf = nCycles * τ

# System solver
systemSolver = create_NLsolver(method=:newton,displayStatus=false)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=pitchingAirfoilDSModelTest,finalTime=tf,Δt=Δt,systemSolver=systemSolver)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
aeroVars = sol.aeroVariables
α = [rad2deg.(aeroVars[i][1].steadyKin.α) for i in eachindex(t)]
cn = [aeroVars[i][1].aeroCoefficients.cn for i in eachindex(t)]
cm = [aeroVars[i][1].aeroCoefficients.cm for i in eachindex(t)]
ct = [aeroVars[i][1].aeroCoefficients.ct for i in eachindex(t)]
cl = @. cn*cosd(α) + ct*sind(α)
cdrag = @. cn*sind(α) - ct*cosd(α)

# Frame string
frameString = source == "NASA" ? 
              string(frame) : 
              string(frame, pad=8)

# Load reference data
clRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/"*source*"frames/"*frameString*"_cl.txt")
cmRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/"*source*"frames/"*frameString*"_cm.txt")
cdRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/"*source*"frames/"*frameString*"_cd.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/pitchingAirfoilDSModelTest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 3
gr()
refLabel = if source == "NASA"
    "Exp - McAlister et al. (1982)"
elseif source == "GU"
    "Exp. - Green & Giuni (2017)"
end

# Range of last cycle
indBeginLastCycle = argmin(abs.(t .- (t[end]-τ)))
rangeLastCycle = indBeginLastCycle:length(t)

# Pitch angle vs time
plt_alpha = plot(xlabel="Time [s]", ylabel="Pitch angle [deg]")
plot!(t, α, color=:black, lw=lw, label=false)
display(plt_alpha)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_alpha.pdf"))

# cn vs time
plt_cnt = plot(xlabel="Time [s]", ylabel="\$c_n\$")
plot!(t, cn, color=:black, lw=lw, label=false)
display(plt_cnt)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cnt.pdf"))

# cm vs time
plt_cmt = plot(xlabel="Time [s]", ylabel="\$c_m\$")
plot!(t, cm, color=:black, lw=lw, label=false)
display(plt_cmt)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cmt.pdf"))

# ct vs time
plt_ctt = plot(xlabel="Time [s]", ylabel="\$c_t\$")
plot!(t, ct, color=:black, lw=lw, label=false)
display(plt_ctt)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_ctt.pdf"))

# cl or cn vs α
if source == "NASA"
    plt_cla = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_l\$")
    plot!(α[rangeLastCycle], cl[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
    scatter!(clRef[1,:], clRef[2,:], color=:black, ms=ms, label=refLabel)
    display(plt_cla)
    savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cla.pdf"))
elseif source == "GU"
    plt_cna = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_n\$")
    plot!(α[rangeLastCycle], cn[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
    scatter!(cnRef[1,:], cnRef[2,:], color=:black, ms=ms, label=refLabel)
    display(plt_cna)
    savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cna.pdf"))
end

# cm vs α
plt_cma = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_m\$", legend=:bottomleft)
plot!(α[rangeLastCycle], cm[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
scatter!(cmRef[1,:], cmRef[2,:], color=:black, ms=ms, label=refLabel)
display(plt_cma)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cma.pdf"))

# cd or ct vs α
if source == "NASA"
    plt_cda = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_d\$")
    plot!(α[rangeLastCycle], cdrag[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
    scatter!(cdRef[1,:], cdRef[2,:], color=:black, ms=ms, label=refLabel)
    display(plt_cda)
    savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cda.pdf"))
elseif source == "GU"
    plt_cta = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_t\$")
    plot!(α[rangeLastCycle], ct[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
    scatter!(ctRef[1,:], ctRef[2,:], color=:black, ms=ms, label=refLabel)
    display(plt_cta)
    savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cta.pdf"))
end

println("Finished pitchingAirfoilDSModelTest.jl")