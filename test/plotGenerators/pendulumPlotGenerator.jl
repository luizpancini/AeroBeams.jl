using Plots

# Run the script
include("../examples/pendulum.jl")

# Print waring
if abs(θ₀) > π/8
    @info "Initial angle of release is large, analytical comparison is not fair"
end

# Set paths
relPath = "/test/outputs/figures/pendulum"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,plotFrequency=1,plotLimits=([-L,L],[-L,L],[-L,0]),save=true,savePath=string(relPath,"/pendulum_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
ms = 4
msw = 0
gr()

# Normalized tip u1 displacement
plt1 = plot(xlabel="\$t/T\$", ylabel="Tip \$u_1/L\$")
plot!(t/T,u1_tip/L, c=:black, lw=lw, label="Numerical")
scatter!(t[1:2:end]/T,u1_tip_analytical[1:2:end]/L, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt1)
savefig(string(absPath,"/pendulum_u1.pdf"))

# Normalized tip u3 displacement
plt2 = plot(xlabel="\$t/T\$", ylabel="Tip \$u_3/L\$")
plot!(t/T,u3_tip/L, c=:black, lw=lw, label="Numerical")
scatter!(t[1:2:end]/T,u3_tip_analytical[1:2:end]/L, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt2)
savefig(string(absPath,"/pendulum_u3.pdf"))

println("Finished pendulumPlotGenerator.jl")