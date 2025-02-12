using Plots, ColorSchemes

# Run the script
include("../examples/PazyFFWTsteadyURangeFixedFold.jl")

# Set paths
relPath = "/test/outputs/figures/PazyFFWTsteadyURangeFixedFold"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape at last airspeed
deformationPlot = plot_steady_deformation(problem[end],plotUndeformed=false,plotDistLoads=false,view=(30,30),legendPos=(0.3,0.5),save=true,savePath=string(relPath,"/PazyFFWTsteadyURangeFixedFold_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Normalized deformed wingspan seen from x2-plane
plt_x1x3 = plot(xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position", xlims=[0,1])
for (i,U) in enumerate(URange)
    plot!((r_n1.+u1[i])/x1_n[end], (r_n3.+u3[i])/x1_n[end], aspect_ratio=:equal, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt_x1x3)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_disp_x2plane.pdf"))

# Normalized deformed wingspan seen from x1-plane
plt_x2x3 = plot(xlabel="Normalized in-plane position", ylabel="Normalized out-of-plane position", xlims=[-0.1, 0.1])
for (i,U) in enumerate(URange)
    plot!((r_n2.+u2[i])/x1_n[end], (r_n3.+u3[i])/x1_n[end], aspect_ratio=:equal, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt_x2x3)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_disp_x1plane.pdf"))

# Normalized deformed wingspan seen from x3-plane
plt_x1x2 = plot(xlabel="Normalized spanwise position", ylabel="Normalized in-plane position", xlims=[0, 1])
for (i,U) in enumerate(URange)
    plot!((r_n1.+u1[i])/x1_n[end], (r_n2.+u2[i])/x1_n[end], aspect_ratio=:equal, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt_x1x2)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_disp_x3plane.pdf"))

# AoA
plt_aoa = plot(xlabel="Normalized spanwise position", ylabel="\$\\alpha\$ [deg]")
for (i,U) in enumerate(URange)
    if U==0
        continue
    end
    plot!(x1_e/x1_n[end], α[i]*180/π, lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt_aoa)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_AoA.pdf"))

# cn
plt_cn = plot(xlabel="Normalized spanwise position", ylabel="\$c_n\$")
for (i,U) in enumerate(URange)
    plot!(x1_e/x1_n[end], cn[i], lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt_cn)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_cn.pdf"))

# Hinge rotation parameters
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
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_deltaphinge.pdf"))

# Hinge angle
plt_phi = plot(xlabel="Airspeed [m/s]", ylabel="\$\\phi\$ [deg]", ylims=[-180,180], yticks=[-180,-90,0,90,180])
plot!(URange, ϕHinge, lw=lw, c=:black, label=false)
display(plt_phi)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_phi.pdf"))

# Hinge moment
plt_M = plot(xlabel="Airspeed [m/s]", ylabel="Hinge moment [N.m]")
plot!(URange, hingeMoment, lw=lw, c=:black, label=false)
display(plt_M)
savefig(string(absPath,"/PazyFFWTsteadyURangeFixedFold_M.pdf"))

println("Finished PazyFFWTsteadyURangeFixedFoldPlotGenerator.jl")