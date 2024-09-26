using Plots, ColorSchemes

# Run the script
include("../examples/conventionalHALEtrim.jl")

# Set paths
relPath = "/test/outputs/figures/conventionalHALEtrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformation plot
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/conventionalHALEtrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
colorsU = get(colorschemes[:rainbow], LinRange(0, 1, length(URange)))
colors = get(colorschemes[:rainbow], LinRange(0, 1, 2))
labels = ["Elastic" "Rigid"]
lw = 2
ms = 3
msw = 0
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root angle of attack [deg]", xlims=[URange[1],URange[end]], ylims=[0,20])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, msw=msw, label="Patil et al. (2001)")
for (i,λ) in enumerate(λRange)
    plot!(URange, trimAoA[i,:], c=colors[i], lw=lw, label=labels[i])
    if i==1
        scatter!(trimAoAERef[1,:], trimAoAERef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
    else
        scatter!(trimAoARRef[1,:], trimAoARRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
    end
end
savefig(string(absPath,"/conventionalHALEtrim_AoA.pdf"))
display(plt1)

# Trim deflected wingspan at U = 25 m/s
U2plot = 25.0
indU = findfirst(x->x==U2plot,URange)
if !isnothing(indU)
    plt2 = plot(xlabel="Spanwise length [m]", ylabel="Vertical displacement [m]", xlims=[0,16], ylims=[0,16], xticks=collect(0:4:16), yticks=collect(0:4:16))
    plot!(x1_0.+(trim_u1[1,indU].-trim_u1[1,indU][1]), x3_0.+trim_u3[1,indU].-trim_u3[1,indU][1], c=:black, lw=lw, label="AeroBeams")
    scatter!(trimDispRef[1,:], trimDispRef[2,:], c=:black, ms=ms, msw=msw, label="Patil et al. (2001)")
    display(plt2)
    savefig(string(absPath,"/conventionalHALEtrim_disp.pdf"))
end

# Trim deflected wingspan over airspeed
plt3 = plot(xlabel="Normalized spanwise length", ylabel="Vertical displacement [% semispan]", xlims=[0,1], ylims=[0,100])
for (i,U) in enumerate(URange)
    plot!((x1_0.+(trim_u1[1,i].-trim_u1[1,i][1]))/16, (x3_0.+(trim_u3[1,i].-trim_u3[1,i][1]))/16*100, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt3)
savefig(string(absPath,"/conventionalHALEtrim_u3OverU.pdf"))

println("Finished conventionalHALEtrimPlotGenerator.jl")