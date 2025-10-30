using Plots, ColorSchemes

# Run the script
include("../examples/clampedHingedBeamSpringRange.jl")

# Set paths
relPath = "/test/outputs/figures/clampedHingedBeamSpringRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, [0, 0.0001, 0.1])
lw = 2
gr()

# Deformed shape
plt_u = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$")
for (i,k) in enumerate(kRange)
    plot!((r_n1.+u1[i])/x1[end], (r_n3.+u3[i])/x1[end], lz=k, c=colors, lw=lw, label=false, colorbar_title="Stiffness [N/m]")
end
display(plt_u)
savefig(string(absPath,"/clampedHingedBeamSpringRange_u.pdf"))

# F3
plt_F3 = plot(xlabel="\$x_1/L\$", ylabel="\$F_3^{\\star}\$ [N]")
for (i,k) in enumerate(kRange)
    plot!(x1/L, F3[i], lz=k, c=colors, lw=lw, label=false, colorbar_title="Stiffness [N/m]")
end
display(plt_F3)
savefig(string(absPath,"/clampedHingedBeamSpringRange_F3.pdf"))

# M2
plt_M2 = plot(xlabel="\$x_1/L\$", ylabel="\$M_2^{\\star}\$ [N.m]")
for (i,k) in enumerate(kRange)
    plot!(x1/L, M2[i], lz=k, c=colors, lw=lw, label=false, colorbar_title="Stiffness [N/m]")
end
display(plt_M2)
savefig(string(absPath,"/clampedHingedBeamSpringRange_M2.pdf"))

# Spring moment
maximumTheoreticalSpringMoment = abs(q₀)*(L/2)*(L/2)/2
plt_Ms = plot(xlabel="\$log(k)\$", ylabel="\$M_s\$ [N.m]", xlims=log10.([minimum(kRange), maximum(kRange)]), ylims=[0,maximumTheoreticalSpringMoment], yticks=0:0.025:maximumTheoreticalSpringMoment)
plot!(log10.(kRange),springMoment, lw=lw, marker=:circle, label=false)
display(plt_Ms)
savefig(string(absPath,"/clampedHingedBeamSpringRange_Ms.pdf"))

# Hinge angle
plt_phi = plot(xlabel="\$log(k)\$", ylabel="Hinge angle [deg]", xlims=log10.([minimum(kRange), maximum(kRange)]), ylims=[0,90], yticks=[-90,-45,0,45,90])
plot!(log10.(kRange),hingeAngle, lw=lw, marker=:circle, label=false)
display(plt_phi)
savefig(string(absPath,"/clampedHingedBeamSpringRange_phi.pdf"))

println("Finished clampedHingedBeamSpringRangePlotGenerator.jl")