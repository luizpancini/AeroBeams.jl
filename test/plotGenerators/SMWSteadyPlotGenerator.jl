using Plots

# Run the script
include("../examples/SMWSteady.jl")

# Set paths
relPath = "/test/outputs/figures/SMWSteady"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/SMWSteady_deformation.pdf"))
display(deformationPlot)

# Plot configurations
gr()
lw = 2
ms = 3

# Normalized deformed wingspan
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$ [% semispan]", xlims=[0,1], ylims=[-20,60])
for (i,U) in enumerate(URange)
    plot!(x1_def[i]/L, x3_def[i]/L*100, lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt1)
savefig(string(absPath,"/SMWSteady_disp.pdf"))

# Tip OOP disp vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp [% semispan]", xlims=[0,30], ylims=[-20,60])
plot!(URange, tip_u3/L*100, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/SMWSteady_OOP.pdf"))

# Tip twist vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]", xlims=[0,30], ylims=[-1,3])
plot!(URange, tip_twist, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/SMWSteady_twist.pdf"))

println("Finished SMWSteady.jl")