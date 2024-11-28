using Plots

# Run the script
include("../examples/PazyFFWTsteadyURange.jl")

# Set paths
relPath = "/test/outputs/figures/PazyFFWTsteadyURange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape at last airspeed
deformationPlot = plot_steady_deformation(problem[end],plotUndeformed=false,plotDistLoads=false,view=(30,30),legendPos=(0.3,0.5),save=true,savePath=string(relPath,"/PazyFFWTsteadyURange_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Normalized deformed wingspan
plt1 = plot(xlabel="Normalized spanwise position", ylabel="Normalized out-of-plane position", xlims=[0,1])
for (i,U) in enumerate(URange)
    plot!((r_n1.+u1[i])/x1_n[end], (r_n3.+u3[i])/x1_n[end], aspect_ratio=:equal, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt1)
savefig(string(absPath,"/PazyFFWTsteadyURange_disp.pdf"))

# AoA
plt2 = plot(xlabel="Normalized spanwise position", ylabel="\$\\alpha\$ [deg]")
for (i,U) in enumerate(URange)
    if U==0
        continue
    end
    plot!(x1_e/x1_n[end], α[i]*180/π, lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt2)
savefig(string(absPath,"/PazyFFWTsteadyURange_AoA.pdf"))

# cn
plt3 = plot(xlabel="Normalized spanwise position", ylabel="\$c_n\$ ")
for (i,U) in enumerate(URange)
    plot!(x1_e/x1_n[end], cn[i], lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt3)
savefig(string(absPath,"/PazyFFWTsteadyURange_cn.pdf"))

println("Finished PazyFFWTsteadyURangePlotGenerator.jl")