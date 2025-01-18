using Plots, ColorSchemes

# Run the script
include("../examples/HealyBaselineFFWTsteadyAoARangeURangeCoast.jl")

# Set paths
relPath = "/test/outputs/figures/HealyBaselineFFWTsteadyAoARangeURangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
lw = 2
ms = 6
msw = 0
gr()

# Hinge OOP displacement
pltu = plot(xlabel="Airspeed [m/s]", ylabel="Hinge displacement [\$\\%\$ of semispan]", xlims=[0,40], ylims=[-10,20], yticks=-10:10:20, legendfontsize=10, legend=:topleft)
plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (i,θ) in enumerate(θRange)
    if i==1
        disp = disp_aoa_25
    elseif i==2
        disp = disp_aoa_50
    else
        disp = disp_aoa_75
    end
    plot!(URange, u3Hinge[i,:]/L*100, lw=lw, ls=:solid, c=colors[i], label="\$\\theta = $(round(θ*180/π,digits=1)) \\degree\$")
    plot!(disp[1,:], disp[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
end
display(pltu)
savefig(string(absPath,"/HealyBaselineFFWTsteadyAoARangeURangeCoast_u3Hinge.pdf"))

# Root bending moment
pltM = plot(xlabel="Airspeed [m/s]", ylabel="Root bending moment [N.m]", xlims=[0,40], ylims=[-10,40], yticks=-10:10:40, legendfontsize=10, legend=:topleft)
plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (i,θ) in enumerate(θRange)
    if i==1
        RBM = RBM_aoa_25
    elseif i==2
        RBM = RBM_aoa_50
    else
        RBM = RBM_aoa_75
    end
    plot!(URange, -M2root[i,:], lw=lw, ls=:solid, c=colors[i], label="\$\\theta = $(round(θ*180/π,digits=1)) \\degree\$")
    plot!(RBM[1,:], RBM[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
end
display(pltM)
savefig(string(absPath,"/HealyBaselineFFWTsteadyAoARangeURangeCoast_M2root.pdf"))

# Lift
pltLift = plot(xlabel="Airspeed [m/s]", ylabel="Lift [N]", xlims=[0,40], ylims=[-10,60], yticks=-10:10:60, legendfontsize=10, legend=:topleft)
plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (i,θ) in enumerate(θRange)
    if i==1
        liftHealy = lift_aoa_25
    elseif i==2
        liftHealy = lift_aoa_50
    else
        liftHealy = lift_aoa_75
    end
    plot!(URange, lift[i,:], lw=lw, ls=:solid, c=colors[i], label="\$\\theta = $(round(θ*180/π,digits=1)) \\degree\$")
    plot!(liftHealy[1,:], liftHealy[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
end
display(pltLift)
savefig(string(absPath,"/HealyBaselineFFWTsteadyAoARangeURangeCoast_lift.pdf"))

# Coast angle
pltϕ = plot(xlabel="Airspeed [m/s]", ylabel="Coast angle [deg]", xlims=[0,40], ylims=[-90,90], yticks=-90:30:90, legendfontsize=10, legend=:bottomright)
plot!([NaN],[NaN], lc=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
plot!([NaN],[NaN], lc=:black, ls=:solid, lw=lw, label="AeroBeams")
for (i,θ) in enumerate(θRange)
    if i==1
        fold = fold_aoa_25
    elseif i==2
        fold = fold_aoa_50
    else
        fold = fold_aoa_75
    end
    plot!(URange, -ϕHinge[i,:], lw=lw, ls=:solid, c=colors[i], label="\$\\theta = $(round(θ*180/π,digits=1)) \\degree\$")
    plot!(fold[1,:], fold[2,:], lw=lw, ls=:dash, c=colors[i], label=false)
end
display(pltϕ)
savefig(string(absPath,"/HealyBaselineFFWTsteadyAoARangeURangeCoast_coastAngle.pdf"))

println("Finished HealyBaselineFFWTsteadyAoARangeURangeCoastPlotGenerator.jl")