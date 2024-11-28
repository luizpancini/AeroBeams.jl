using Plots 

# Run the script
include("../examples/clampedHingedBeamSpringRange.jl")

# Set paths
relPath = "/test/outputs/figures/clampedHingedBeamSpringRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 1
gr()

# u3
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
for (i,k) in enumerate(kRange)
    plot!((r_n1.+u1[i])/x1[end], (r_n3.+u3[i])/x1[end], lz=k, c=:rainbow, lw=lw, label=false, colorbar_title="Stiffness [N/m]")
end
display(plt1)
savefig(string(absPath,"/clampedHingedBeamSpringRange_u3.pdf"))

# p2
plt2 = plot(xlabel="\$x_1/L\$", ylabel="\$p_2\$")
for (i,k) in enumerate(kRange)
    plot!(x1/L, p2[i], lz=k, c=:rainbow, lw=lw, label=false, colorbar_title="Stiffness [N/m]")
end
display(plt2)
savefig(string(absPath,"/clampedHingedBeamSpringRange_p2.pdf"))

# F3
plt3 = plot(xlabel="\$x_1/L\$", ylabel="\$F_3\$ [N]")
for (i,k) in enumerate(kRange)
    plot!(x1/L, F3[i], lz=k, c=:rainbow, lw=lw, label=false, colorbar_title="Stiffness [N/m]")
end
display(plt3)
savefig(string(absPath,"/clampedHingedBeamSpringRange_F3.pdf"))

# M2
plt4 = plot(xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]")
for (i,k) in enumerate(kRange)
    plot!(x1/L, M2[i], lz=k, c=:rainbow, lw=lw, label=false, colorbar_title="Stiffness [N/m]")
end
display(plt4)
savefig(string(absPath,"/clampedHingedBeamSpringRange_M2.pdf"))

println("Finished clampedHingedBeamSpringRangePlotGenerator.jl")