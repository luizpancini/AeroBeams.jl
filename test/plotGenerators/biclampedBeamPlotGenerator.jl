using Plots

# Run the script
include("../examples/biclampedBeam.jl")

# Set paths
relPath = "/test/outputs/figures/biclampedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,scale=1/δ,plotLimits=([0,L],[-0.25,0.25],[0,1]),save=true,savePath=string(relPath,"/biclampedBeam_deformation.gif"))

# Plot configurations
lw = 2
ms = 6
msw = 0
gr()

# Normalized displacement at mid-length
plt1 = plot(xlabel="\$t/T\$", ylabel="\$u_3/\\delta\$ at \$x_1=L/2\$")
plot!(tNorm,u3_mid/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],u3_mid_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt1)
savefig(string(absPath,"/biclampedBeam_midDisp.pdf"))

# Normalized velocity at mid-length 
plt2 = plot(xlabel="\$t/T\$", ylabel="\$V_3/\\delta\$ at \$x_1=L/2\$ [\$1\$/s]")
plot!(tNorm,V3_mid/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],V3_mid_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt2)
savefig(string(absPath,"/biclampedBeam_midVel.pdf"))

# Normalized acceleration at mid-length
plt3 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3/\\delta\$ at \$x_1=L/2\$ [\$1\$/\$s^2\$]")
plot!(tNorm,Vdot3_mid/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Vdot3_mid_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt3)
savefig(string(absPath,"/biclampedBeam_midAcc.pdf"))

# Normalized rotation at quarter-length
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\theta/(2\\pi\\delta)\$ at \$x_1=L/4\$")
plot!(tNorm,θ2_quarter/(2*π)/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],θ2_quarter_analytic[1:5:end]/(2*π)/δ, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt4)
savefig(string(absPath,"/biclampedBeam_quarterRot.pdf"))

# Angular velocity at quarter-length
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_2\$ at \$x_1=L/4\$ [rad/s]")
plot!(tNorm,Ω2_quarter, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ω2_quarter_analytic[1:5:end], c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt5)
savefig(string(absPath,"/biclampedBeam_quarterAngVel.pdf"))

# Angular acceleration at quarter-length
plt6 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_2\$ at \$x_1=L/4\$ [rad/\$s^2\$]")
plot!(tNorm,Ωdot2_quarter, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ωdot2_quarter_analytic[1:5:end], c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt6)
savefig(string(absPath,"/biclampedBeam_quarterAngAcc.pdf"))

println("Finished biclampedBeamPlotGenerator.jl")