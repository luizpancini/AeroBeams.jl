using Plots

# Run the script
include("../examples/initialDisplacementBeam.jl")

# Set paths
relPath = "/test/outputs/figures/initialDisplacementBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,scale=1/δ/10,scalePos=[0.15;-0.05;0],timeStampPos=[0.5;0.05;0],plotFrequency=1,plotLimits=([0,L],[-L/2,L/2],[-L/10,L/10]),save=true,savePath=string(relPath,"/initialDisplacementBeam_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
ms = 5
gr()

# Normalized displacement at quarter-length
plt1 = plot(xlabel="\$t/T\$", ylabel="\$u_3/\\delta\$ at \$x_1=L/4\$")
plot!(tNorm,u3_quarter/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],u3_quarter_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt1)
savefig(string(absPath,"/initialDisplacementBeam_disp.pdf"))

# Normalized velocity at quarter-length
plt2 = plot(xlabel="\$t/T\$", ylabel="\$V_3/\\delta\$ at \$x_1=L/4\$ [\$1\$/s]")
plot!(tNorm,V3_quarter/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],V3_quarter_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt2)
savefig(string(absPath,"/initialDisplacementBeam_vel.pdf"))

# Normalized acceleration at quarter-length
plt3 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3/\\delta\$ at \$x_1=L/4\$ [\$1\$/\$s^2\$]")
plot!(tNorm,Vdot3_quarter/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Vdot3_quarter_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt3)
savefig(string(absPath,"/initialDisplacementBeam_acc.pdf"))

# Normalized rotation at root
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\theta/(2\\pi\\delta)\$ at \$x_1=0\$")
plot!(tNorm,θ2_root/(2*π)/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],θ2_root_analytic[1:5:end]/(2*π)/δ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt4)
savefig(string(absPath,"/initialDisplacementBeam_rot.pdf"))

# Angular velocity at mid-length
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_2\$ at \$x_1=L/2\$ [rad/s]")
plot!(tNorm,Ω2_mid, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ω2_mid_analytic[1:5:end], c=:blue, ms=ms, msw=0, label="Analytical")
display(plt5)
savefig(string(absPath,"/initialDisplacementBeam_angVel.pdf"))

# Angular acceleration at mid-length
plt6 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_2\$ at \$x_1=L/2\$ [rad/\$s^2\$]")
plot!(tNorm,Ωdot2_mid, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ωdot2_mid_analytic[1:5:end], c=:blue, ms=ms, msw=0, label="Analytical")
display(plt6)
savefig(string(absPath,"/initialDisplacementBeam_angAcc.pdf"))

println("Finished initialDisplacementBeamPlotGenerator.jl")