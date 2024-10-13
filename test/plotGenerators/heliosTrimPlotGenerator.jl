using Plots, ColorSchemes

# Run the script
include("../examples/heliosTrim.jl")

# Set paths
relPath = "/test/outputs/figures/heliosTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape of flexible aircraft at highest payload
deformationPlot = plot_steady_deformation(problem[1,end],view=(30,30),save=true,savePath=string(relPath,"/heliosTrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(λRange)))
labels = ["Flexible" "Rigid"]
lw = 2
ms = 3
gr()

# Trim root angle of attack
plt1 = plot(xlabel="Payload [lb]", ylabel="Trim root AoA [deg]", xlims=[0,500], ylims=[0,5])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, label="Patil & Hodges (2006)")
for (i,λ) in enumerate(λRange)
    plot!(PRange, trimAoA[i,:], c=colors[i], lw=lw, label=labels[i])
    if i==1
        scatter!(αFlexibleRef[1,:], αFlexibleRef[2,:], c=colors[i], ms=ms, msw=0, label=false)
    else
        scatter!(αRigidRef[1,:], αRigidRef[2,:], c=colors[i], ms=ms, msw=0, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/heliosTrim_AoA.pdf"))

# Trim thrust force
plt2 = plot(xlabel="Payload [lb]", ylabel="Trim thrust per motor [N]", xlims=[0,500], ylims=[0,40])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, label="Patil & Hodges (2006)")
for (i,λ) in enumerate(λRange)
    plot!(PRange, trimThrust[i,:], c=colors[i], lw=lw, label=labels[i])
    scatter!(TRef[1,:], TRef[2,:], c=colors[i], ms=ms, msw=0, label=false)
end
display(plt2)
savefig(string(absPath,"/heliosTrim_thrust.pdf"))

# Trim elevator deflection
plt3 = plot(xlabel="Payload [lb]", ylabel="Trim elevator deflection [deg]", xlims=[0,500], ylims=[0,10])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, label="Patil & Hodges (2006)")
for (i,λ) in enumerate(λRange)
    plot!(PRange, trimδ[i,:], c=colors[i], lw=lw, label=labels[i])
    if i==1
        scatter!(δFlexibleRef[1,:], δFlexibleRef[2,:], c=colors[i], ms=ms, msw=0, label=false)
    else
        scatter!(δRigidRef[1,:], δRigidRef[2,:], c=colors[i], ms=ms, msw=0, label=false)
    end
end
display(plt3)
savefig(string(absPath,"/heliosTrim_delta.pdf"))

println("Finished heliosTrimPlotGenerator.jl")