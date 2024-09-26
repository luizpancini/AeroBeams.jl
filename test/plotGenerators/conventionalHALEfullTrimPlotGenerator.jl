using Plots

# Run the script
include("../examples/conventionalHALEfullTrim.jl")

# Set paths
relPath = "/test/outputs/figures/conventionalHALEfullTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformation plot
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/conventionalHALEtrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[20,35], ylims=[0,15])
plot!(URange, trimAoA, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/conventionalHALEfullTrim_AoA.pdf"))

# Trim propeller force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[20,35])
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/conventionalHALEfullTrim_thrust.pdf"))

# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[20,35])
plot!(URange, trimδ, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/conventionalHALEfullTrim_delta.pdf"))

println("Finished conventionalHALEfullTrim.jl")