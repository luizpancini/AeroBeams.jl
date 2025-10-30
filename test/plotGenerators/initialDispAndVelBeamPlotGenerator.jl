using Plots

# Run the script
include("../examples/initialDispAndVelBeam.jl")

# Set paths
relPath = "/test/outputs/figures/initialDispAndVelBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,scale=1/δ/10,scalePos=[0.15;-0.05;0],timeStampPos=[0.5;-0.05;0],plotFrequency=1,plotLimits=([0,L],[-L/2,L/2],[-L/3,L/3]),save=true,savePath=string(relPath,"/initialDispAndVelBeam_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
ts = 10
fs = 16
lw = 2
ms = 5
gr()

# Displacement at quarter-length
plt1 = plot(xlabel="\$t/T\$", ylabel="\$u_3\$ at \$x_1=L/4\$ [m]", tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
plot!(tNorm,u3_quarter, c=:black, lw=lw, label="AeroBeams")
scatter!(tNorm[1:5:end],u3_quarter_analytic[1:5:end], c=:black, ms=ms, msw=0, label="Analytical")
display(plt1)
savefig(string(absPath,"/initialDispAndVelBeam_disp.pdf"))

# Velocity at quarter-length
plt2 = plot(xlabel="\$t/T\$", ylabel="\$V_3^\\star\$ at \$x_1=L/4\$ [m/s]", tickfont=font(ts), guidefont=font(fs))
plot!(tNorm,V3_quarter, c=:black, lw=lw, label=false)
scatter!(tNorm[1:5:end],V3_quarter_analytic[1:5:end], c=:black, ms=ms, msw=0, label=false)
display(plt2)
savefig(string(absPath,"/initialDispAndVelBeam_vel.pdf"))

# Acceleration at quarter-length
plt3 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3^\\star\$ at \$x_1=L/4\$ [m/\$s^2\$]", tickfont=font(ts), guidefont=font(fs))
plot!(tNorm,Vdot3_quarter, c=:black, lw=lw, label=false)
scatter!(tNorm[1:5:end],Vdot3_quarter_analytic[1:5:end], c=:black, ms=ms, msw=0, label=false)
display(plt3)
savefig(string(absPath,"/initialDispAndVelBeam_acc.pdf"))

# Rotation at root
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\theta/(2\\pi)\$ at \$x_1=0\$ [rad]", tickfont=font(ts), guidefont=font(fs))
plot!(tNorm,θ2_root/(2*π), c=:black, lw=lw, label=false)
scatter!(tNorm[1:5:end],θ2_root_analytic[1:5:end]/(2*π), c=:black, ms=ms, msw=0, label=false)
display(plt4)
savefig(string(absPath,"/initialDispAndVelBeam_rot.pdf"))

# Angular velocity at mid-length
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_2^\\star\$ at \$x_1=L/2\$ [rad/s]", tickfont=font(ts), guidefont=font(fs))
plot!(tNorm,Ω2_mid, c=:black, lw=lw, label=false)
scatter!(tNorm[1:5:end],Ω2_mid_analytic[1:5:end], c=:black, ms=ms, msw=0, label=false)
display(plt5)
savefig(string(absPath,"/initialDispAndVelBeam_angVel.pdf"))

# Angular acceleration at mid-length
plt6 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_2^\\star\$ at \$x_1=L/2\$ [rad/\$s^2\$]", tickfont=font(ts), guidefont=font(fs))
plot!(tNorm,Ωdot2_mid, c=:black, lw=lw, label=false)
scatter!(tNorm[1:5:end],Ωdot2_mid_analytic[1:5:end], c=:black, ms=ms, msw=0, label=false)
display(plt6)
savefig(string(absPath,"/initialDispAndVelBeam_angAcc.pdf"))

# Acceleration of elements
plot_time_outputs(problem,elements=collect(1:nElem),elementalOutputs=["Vdot3"],save=true,saveFolder=string(relPath,"/"))

println("Finished initialDispAndVelBeamPlotGenerator.jl")