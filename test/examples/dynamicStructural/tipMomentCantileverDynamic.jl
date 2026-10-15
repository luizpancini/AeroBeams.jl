using AeroBeams

# Beam
L = EIy = 1
ρA = ρI = 1e-6
nElem = 40
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)],I=[inertia_matrix(ρA=ρA,ρIy=ρI)])

# Time variables
tf = 1
Nt = 100
Δt = tf/Nt

# BCs
M₀ = 2π*EIy/L
M = t -> M₀*t/tf
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipMoment = create_BC(name="tipMoment",beam=beam,node=nElem+1,conditions=Dict(:M2A=>t->M(t)))

# Model
tipMomentCantileverDynamic = create_Model(name="tipMomentCantileverDynamic",beams=[beam],BCs=[clamp,tipMoment],units=create_UnitsSystem(length="in",force="lbf"))

# Create and solve the problem
problem = create_DynamicProblem(model=tipMomentCantileverDynamic,finalTime=tf,Δt=Δt)
solve!(problem)

# Get solution over time
sol = problem.solutionHistory
t = sol.t
states = sol.states
κ2 = [states[i][1].midpoint.κ[2] for i in eachindex(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/tipMomentCantileverDynamic"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ts = 20
fs = 24
lw = 2
fps = 5
pltFreq = round(Int,tf/Δt/20)
DPI = 300

# Animation
anim = plot_dynamic_deformation(problem,plotAxes=false,plotFrequency=pltFreq,DPI=DPI,lw=lw,fps=fps,plotLimits=([-L/2,L],[0,L],[-3L/4,3L/4]),view=(0,0),save=true,savePath=string(relPath,"/tipMomentCantileverDynamic_deformation.gif"),displayProgress=true)
display(anim)

# Animated plot of curvature
gr()
tSpaced = t[1:pltFreq:end]
κ2Spaced = κ2[1:pltFreq:end]
plt = plot(xlabel="\$M_2\$", ylabel="\$\\kappa_2\$", title="\$\\kappa_2 = K_2 - k_2\$", xlims=[0,tf], ylims=[0, 2π], yticks=([0, 2π], ["\$0\$", "\$2π\$"]), xticks=false, tickfont=font(ts), guidefont=font(fs), titlefont=font(fs), dpi=DPI)
anim = @animate for (tt,timeNow) in enumerate(tSpaced)
    plot!(tSpaced[1:tt], κ2Spaced[1:tt], c=:black, label=false)
end
gif_handle = gif(anim, string(absPath,"/tipMomentCantileverDynamicPlotGenerator_kappa2.gif"), fps=fps)
display(gif_handle)

println("Finished tipMomentCantileverDynamic.jl")
