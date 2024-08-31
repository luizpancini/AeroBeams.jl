using Plots

# Run the script
include("../examples/flapOscillationLinked.jl")

# Set paths
relPath = "/test/outputs/figures/flapOscillationLinked"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=10,plotLimits=[(0,2*L),(-L,L),(-L,L)],save=true,savePath=string(relPath,"/flapOscillationLinked_deformation.gif"),displayProgress=true,plotAeroSurf=false)

# Plot configurations
lw = 2
labels = ["Master" "Slave"]
gr()

# cn and cm vs time
plt11 = plot(ylabel="\$c_n\$", xlims=[0,cycles])
plot!(tNorm, [cnMaster, cnSlave], lw=lw, label=labels)
plt12 = plot(xlabel="\$t/T\$", ylabel="\$c_m\$", xlims=[0,cycles])
plot!(tNorm, [cmMaster, cmSlave], lw=lw, label=false)
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/flapOscillationLinked.pdf"))

println("Finished flapOscillationLinkedPlotGenerator.jl")