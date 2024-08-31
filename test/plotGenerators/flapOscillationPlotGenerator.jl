using Plots

# Run the script
include("../examples/flapOscillation.jl")

# Set paths
relPath = "/test/outputs/figures/flapOscillation"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=10,plotLimits=[(0,L),(-L/2,L/2),(-L/2,L/2)],save=true,savePath=string(relPath,"/flapOscillation_deformation.gif"),displayProgress=true,plotAeroSurf=false)

# Plot configurations
rangeLastCycle = findfirst(x->x==tf-T,t):length(t)
lw = 2
ms = 3
gr()

# cn vs δ
plt1 = plot(xlabel="\$\\delta\$ [deg]", ylabel="\$c_n/\\pi\$", xlims=[-3,3], ylims=[-0.075,0.075])
plot!(δ.(t[rangeLastCycle])*180/π, cn[rangeLastCycle]/π, c=:black, lw=lw, label="AeroBeams")
plot!(cnRefMod[1,:], cnRefMod[2,:]/π, c=:black, ls=:dash, lw=lw, label="Incompressible thoery by Leishman (2006)")
scatter!(cnExp[1,:], cnExp[2,:]/π, c=:black, ms=ms, msw=0, label="Experiment by Tijdeman & Schippers (1973)")
display(plt1)
savefig(string(absPath,"/flapOscillation_cn.pdf"))

# cm vs δ
plt2 = plot(xlabel="\$\\delta\$ [deg]", ylabel="\$-2c_m/\\pi\$", xlims=[-3,3], ylims=[-0.03,0.03])
plot!(δ.(t[rangeLastCycle])*180/π, -2*cm[rangeLastCycle]/π, c=:black, lw=lw, label="AeroBeams")
plot!(cmRefMod[1,:], -2*cmRefMod[2,:]/π, c=:black, ls=:dash, lw=lw, label="Incompressible thoery by Leishman (2006)")
scatter!(cmExp[1,:], -2*cmExp[2,:]/π, c=:black, ms=ms, msw=0, label="Experiment by Tijdeman & Schippers (1973)")
display(plt2)
savefig(string(absPath,"/flapOscillation_cm.pdf"))

println("Finished flapOscillation.jl")