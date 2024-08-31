using Plots

# Run the script
include("../examples/tipLoadedCantileverTrim.jl")

# Set paths
relPath = "/test/outputs/figures/tipLoadedCantileverTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 4
msw = 0
gr()

# u3
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$u_3 / (FL/3EI)\$")
plot!(x1/L, u3/(F*L/(3*EI)), c=:black, lw=lw, label="Numerical")
scatter!(x1/L, u3_analytical(x1)/(F*L/(3*EI)), c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt1)
savefig(string(absPath,"/tipLoadedCantileverTrim_u3.pdf"))

# F3
plt2 = plot(xlabel="\$x_1/L\$", ylabel="\$F_3/F\$ [N]", ylims=[-2,0])
plot!(x1/L, F3/F, c=:black, lw=lw, label="Numerical")
scatter!(x1/L, F3_analytical.(x1)/F, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt2)
savefig(string(absPath,"/tipLoadedCantileverTrim_F3.pdf"))

# M2
plt3 = plot(xlabel="\$x_1/L\$", ylabel="\$M_2/(FL)\$")
plot!(x1/L, M2/(F*L), c=:black, lw=2, label="Numerical")
scatter!(x1/L, M2_analytical(x1)/(F*L), c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt3)
savefig(string(absPath,"/tipLoadedCantileverTrim_M2.pdf"))

println("Finished tipLoadedCantileverTrimPlotGenerator.jl")