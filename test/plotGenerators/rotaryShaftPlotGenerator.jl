using Plots

# Run the script
include("../examples/rotaryShaft.jl")

# Set paths
relPath = "/test/outputs/figures/rotaryShaft"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 3
gr()

# Normalized rotation parameter
plt1 = plot(xlabel="\$t/T\$", ylabel="\$p_1/\\Delta\\theta\$ ")
plot!(tNorm,pNum/Δθ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],p.(t[1:20:end])/Δθ, c=:blue, ms=ms, label="Analytical")
display(plt1)
savefig(string(absPath,"/rotaryShaft_p.pdf"))

# Normalized rotation parameter rate 
plt2 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{p}_1/\\Delta\\theta\\omega\$ [\$1\$/s]")
plot!(tNorm,pdotNum/(Δθ*ω), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],pdot.(t[1:20:end])/(Δθ*ω), c=:blue, ms=ms, label="Analytical")
display(plt2)
savefig(string(absPath,"/rotaryShaft_pdot.pdf"))

# Normalized sectional angular velocity 
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_1/\\Delta\\theta\\omega\$ [1/s]")
plot!(tNorm,ΩNum/(Δθ*ω), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],θdot.(t[1:20:end])/(Δθ*ω), c=:blue, ms=ms, label="Analytical")
display(plt4)
savefig(string(absPath,"/rotaryShaft_angVel.pdf"))

# Normalized sectional angular acceleration 
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_1/\\Delta\\theta\\omega^2\$ [1/\$s^2\$]")
plot!(tNorm,ΩdotNum/(Δθ*ω^2), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],θddot.(t[1:20:end])/(Δθ*ω^2), c=:blue, ms=ms, label="Analytical")
display(plt5)
savefig(string(absPath,"/rotaryShaft_angAcc.pdf"))

# Normalized driving torque
plt6 = plot(xlabel="\$t/T\$", ylabel="\$M_1^*\$ [N.m]")
plot!(tNorm,MNum, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:20:end],-θddot.(t[1:20:end])*(ρ*Is), c=:blue, ms=ms, label="Analytical")
display(plt6)
savefig(string(absPath,"/rotaryShaft_torque.pdf"))

println("Finished rotaryShaftPlotGenerator.jl")