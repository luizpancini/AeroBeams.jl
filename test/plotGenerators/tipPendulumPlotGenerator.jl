using Plots

# Run the script
include("../examples/tipPendulum.jl")

# Print warning, if applicable
if abs(θ₀) > π/8
    println("Initial angle of release is large, analytical comparison is not valid")
end

# Set paths
relPath = "/test/outputs/figures/tipPendulum"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,plotDistLoads=false,plotLimits=([-L,L],[-L,L],[-L,0]),save=true,savePath=string(relPath,"/tipPendulum_deformation.gif"),displayProgress=true)

# Plot configurations
gr()
lw = 2
ms = 4
msw = 0

# Normalized tip u1 displacement
plt1 = plot(xlabel="\$t/T\$", ylabel="Tip \$u_1/L\$ ")
plot!(t/T,u1_tip/L, c=:black, lw=lw, label="Numerical")
scatter!(t[1:2:end]/T,u1_tip_analytical[1:2:end]/L, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt1)
savefig(string(absPath,"/tipPendulum_u1.pdf"))

# Normalized tip u3 displacement
plt2 = plot(xlabel="\$t/T\$", ylabel="Tip \$u_3/L\$ ")
plot!(t/T,u3_tip/L, c=:black, lw=lw, label="Numerical")
scatter!(t[1:2:end]/T,u3_tip_analytical[1:2:end]/L, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt2)
savefig(string(absPath,"/tipPendulum_u3.pdf"))

println("Finished tipPendulumPlotGenerator.jl")