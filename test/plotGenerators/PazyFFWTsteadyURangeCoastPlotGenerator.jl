using Plots, ColorSchemes

# Run the script
include("../examples/PazyFFWTsteadyURangeCoast.jl")

# Set paths
relPath = "/test/outputs/figures/PazyFFWTsteadyURangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape at last airspeed
deformationPlot = plot_steady_deformation(problem[end],plotUndeformed=false,plotDistLoads=false,view=(30,30),legendPos=(0.3,0.5),save=true,savePath=string(relPath,"/PazyFFWTsteadyURangeCoast_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Normalized deformed wingspan seen from x2-plane
plt1 = plot(xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position", xlims=[0,1])
for (i,U) in enumerate(URange)
    plot!((r_n1.+u1[i])/x1_n[end], (r_n3.+u3[i])/x1_n[end], aspect_ratio=:equal, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt1)
savefig(string(absPath,"/PazyFFWTsteadyURangeCoast_disp_x2plane.pdf"))

# Normalized deformed wingspan seen from x1-plane
plt11 = plot(xlabel="Normalized in-plane position", ylabel="Normalized out-of-plane position", xlims=[-0.1, 0.1])
for (i,U) in enumerate(URange)
    plot!((r_n2.+u2[i])/x1_n[end], (r_n3.+u3[i])/x1_n[end], aspect_ratio=:equal, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt11)
savefig(string(absPath,"/PazyFFWTsteadyURangeCoast_disp_x1plane.pdf"))

# Normalized deformed wingspan seen from x3-plane
plt12 = plot(xlabel="Normalized spanwise position", ylabel="Normalized in-plane position", xlims=[0, 1])
for (i,U) in enumerate(URange)
    plot!((r_n1.+u1[i])/x1_n[end], (r_n2.+u2[i])/x1_n[end], aspect_ratio=:equal, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt12)
savefig(string(absPath,"/PazyFFWTsteadyURangeCoast_disp_x3plane.pdf"))

# AoA
plt2 = plot(xlabel="Normalized spanwise position", ylabel="\$\\alpha\$ [deg]")
for (i,U) in enumerate(URange)
    if U==0
        continue
    end
    plot!(x1_e/x1_n[end], α[i]*180/π, lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt2)
savefig(string(absPath,"/PazyFFWTsteadyURangeCoast_AoA.pdf"))

# cn
plt3 = plot(xlabel="Normalized spanwise position", ylabel="\$c_n\$")
for (i,U) in enumerate(URange)
    plot!(x1_e/x1_n[end], cn[i], lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt3)
savefig(string(absPath,"/PazyFFWTsteadyURangeCoast_cn.pdf"))

# pHinge
labels = ["\$\\Delta p_1\$" "\$\\Delta p_2\$" "\$\\Delta p_3\$"]
colors = get(colorschemes[:rainbow], LinRange(0, 1, 3))
p1Hinge = [pHinge[i][1] for i in eachindex(URange)]
p2Hinge = [pHinge[i][2] for i in eachindex(URange)]
p3Hinge = [pHinge[i][3] for i in eachindex(URange)]
plt4 = plot(xlabel="Airspeed [m/s]", ylabel="\$\\Delta p\$")
plot!(URange, p1Hinge, lw=lw, c=colors[1], label=labels[1])
plot!(URange, p2Hinge, lw=lw, c=colors[2], label=labels[2])
plot!(URange, p3Hinge, lw=lw, c=colors[3], label=labels[3])
display(plt4)
savefig(string(absPath,"/PazyFFWTsteadyURangeCoast_deltaphinge.pdf"))

# Coast angle
plt5 = plot(xlabel="Airspeed [m/s]", ylabel="Coast angle [deg]", ylims=[-90,90], yticks=-90:45:90)
plot!(URange, -ϕHinge, lw=lw, c=:black, label=false)
display(plt5)
savefig(string(absPath,"/PazyFFWTsteadyURangeCoast_deltaphihinge.pdf"))

println("Finished PazyFFWTsteadyURangeCoastPlotGenerator.jl")