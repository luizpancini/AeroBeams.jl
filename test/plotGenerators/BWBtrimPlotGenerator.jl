using Plots

# Run the script
include("../examples/BWBtrim.jl")

# Set paths
relPath = "/test/outputs/figures/BWBtrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,ΔuDef=[0;0;-167],view=(45,30),save=true,savePath=string(relPath,"/BWBtrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
ms = 4
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[URange[1],URange[end]])
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="UM/NAST")
plot!(URange, trimAoA, c=:black, lw=lw, label=false)
scatter!(trimAoARef[1,:],trimAoARef[2,:], mc=:black, ms=ms, label=false)
display(plt1)
savefig(string(absPath,"/BWBtrim_AoA.pdf"))

# Trim propeller force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[URange[1],URange[end]])
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="UM/NAST")
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
scatter!(trimThrustRef[1,:],trimThrustRef[2,:], mc=:black, ms=ms, label=false)
display(plt2)
savefig(string(absPath,"/BWBtrim_thrust.pdf"))

# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[URange[1],URange[end]])
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="UM/NAST")
plot!(URange, trimδ, c=:black, lw=lw, label=false)
scatter!(trimδRef[1,:],trimδRef[2,:], mc=:black, ms=ms, label=false)
display(plt3)
savefig(string(absPath,"/BWBtrim_delta.pdf"))

println("Finished BWBtrimPlotGenerator.jl")