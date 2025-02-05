using Plots

# Run the script
include("../examples/flyingScissors.jl")

# Set paths
relPath = "/test/outputs/figures/flyingScissors"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="I",plotFrequency=2,plotLimits=([0,2*L],[-L/2,L/2],[-L,0]),save=true,savePath=string(relPath,"/flyingScissors_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
labels = ["Tip A" "Hinge" "Tip B"]
gr()

# Nomalized u1 displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="\$u_1/L\$")
plot!(t,[u1_tipA/L, u1_hinge/L, u1_tipB/L], lw=lw, label=labels)
display(plt1)
savefig(string(absPath,"/flyingScissors_u1.pdf"))

# Nomalized u3 displacements
plt2 = plot(xlabel="\$t\$ [s]", ylabel="\$u_3/L\$")
plot!(t,[u3_tipA/L, u3_hinge/L, u3_tipB/L], lw=lw, label=labels)
display(plt2)
savefig(string(absPath,"/flyingScissors_u3.pdf"))

println("Finished flyingScissorsPlotGenerator.jl")