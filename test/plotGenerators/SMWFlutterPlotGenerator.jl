using Plots, ColorSchemes

# Run the script
include("../examples/SMWFlutter.jl")

# Set paths
relPath = "/test/outputs/figures/SMWFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes
modesPlot = plot_mode_shapes(problem[end],scale=5,view=(30,30),legendPos=(0.3,0.6),save=true,savePath=string(relPath,"/SMWFlutter_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
gr()

# Normalized deformed wingspan
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$", xlims=[0,1])
for (i,U) in enumerate(URange)
    plot!(x1_def[i]/L, x3_def[i]/L, lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt1)
savefig(string(absPath,"/SMWFlutter_disp.pdf"))

# Angle of attack over wingspan
plt11 = plot(xlabel="\$x_1/L\$", ylabel="\$\\alpha\$ [deg]", xlims=[0,1])
for (i,U) in enumerate(URange)
    plot!(x1_e/L, α_of_x1[i]*180/pi, lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt11)
savefig(string(absPath,"/SMWFlutter_AoA.pdf"))

# Root locus
plt2 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-5,1], ylims=[0,50])
for mode in 1:nModes
    plot!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], lw=lw, label="Mode $mode")
end
display(plt2)
savefig(string(absPath,"/SMWFlutter_rootlocus.pdf"))

# V-g-f
plt31 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode], c=modeColors[mode], lw=lw,  label=false)
end
plt32 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.15,0.05])
for mode in 1:nModes
    plot!(URange, modeDampingRatios[mode], c=modeColors[mode], lw=lw,  label="Mode $mode")
end
plt3 = plot(plt31,plt32, layout=(2,1))
display(plt3)
savefig(string(absPath,"/SMWFlutter_Vgf.pdf"))

# Mode shapes at specified velocity
U2plot = 25
ind = findfirst(x->x==U2plot,URange)
plt4 = plot(xlabel="\$x_1/L\$", ylabel="\$u_2\$", title="Chordwise bending mode shape at U=$(URange[ind]) m/s")
plt5 = plot(xlabel="\$x_1/L\$", ylabel="\$u_3\$", title="Flapwise bending mode shape at U=$(URange[ind]) m/s")
plt6 = plot(xlabel="\$x_1/L\$", ylabel="\$p_1\$", title="Torsional mode shape at U=$(URange[ind]) m/s")
for m in 1:nModes
    plot!(plt4,x1_e/L, u2_modeShapes[ind,m], c=modeColors[m], lw=lw, label="Mode $m")
    plot!(plt5,x1_e/L, u3_modeShapes[ind,m], c=modeColors[m], lw=lw, label="Mode $m")
    plot!(plt6,x1_e/L, p1_modeShapes[ind,m], c=modeColors[m], lw=lw, label="Mode $m")
end
display(plt4)
display(plt5)
display(plt6)
savefig(string(absPath,"/SMWFlutter_modalDisp.pdf"))

println("Finished SMWFlutter.jl")