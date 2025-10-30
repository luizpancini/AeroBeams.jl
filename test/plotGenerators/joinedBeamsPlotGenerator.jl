using Plots

# Run the script
include("../examples/joinedBeams.jl")

# Set paths
relPath = "/test/outputs/figures/joinedBeams"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,scale=100,plotFrequency=10,view=(90,0),plotLimits=([-L,L],[-15,15],[0,8]),save=true,savePath=string(relPath,"/joinedBeams_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Plot u1 of left tip over time
plt1 = plot(t, u1, lw=lw, label=false, xlabel="\$t\$ [s]", ylabel="\$u_1\$ [m]")
display(plt1)
savefig(string(absPath,"/joinedBeams_u1.pdf"))

# Plot u2 of left tip over time
plt2 = plot(t, u2, lw=lw, label=false, xlabel="\$t\$ [s]", ylabel="\$u_2\$ [m]")
display(plt2)
savefig(string(absPath,"/joinedBeams_u2.pdf"))

# Plot u3 of left tip over time
plt3 = plot(t, u3, lw=lw, label=false, xlabel="\$t\$ [s]", ylabel="\$u_3\$ [m]")
display(plt3)
savefig(string(absPath,"/joinedBeams_u3.pdf"))

println("Finished joinedBeamsPlotGenerator.jl")