using Plots, ColorSchemes

# Run the script
include("../examples/PazyWingPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed state of last problem
deformationPlot = plot_steady_deformation(problem,view=(45,30),save=true,savePath="/test/outputs/figures/PazyWingPitchRange/PazyWingPitchRange_deformation.pdf")
display(deformationPlot)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(θRange)))
ts = 10
fs = 16
lfs = 9
lw = 2
ms = 3
msw = 0
gr()

# Tip OOP displacement vs. airspeed for root several pitch angles 
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,60], ylims=[0,50], tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
plot!([NaN], [NaN], c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dashdot, label="UM/NAST (Exp. loss) - Riso & Cesnik (2022,2023)")
scatter!([NaN], [NaN], c=:black, ms=ms, label="Exp. (test 1) - Avin et al. (2022)")
scatter!([NaN], [NaN], c=:black, shape=:diamond, ms=ms, label="Exp. (test 2) - Avin et al. (2022)")
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_OOP[i,:]/L*100, c=colors[i], lw=lw, ls=:solid, label=string("\$\\alpha_r=",round(Int,θ),"^\\circ\$"))
    if θ==3
        scatter!(sqrt.(2*tip_u3Vsq_rootPitch3_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch3_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    elseif θ==5
        plot!(tip_u3VsU_rootPitch5_UMNAST[1,:], tip_u3VsU_rootPitch5_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_u3VsU_rootPitch5_Exp[1,:], tip_u3VsU_rootPitch5_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        scatter!(sqrt.(2*tip_u3Vsq_rootPitch5_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch5_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    elseif θ==7
        plot!(tip_u3VsU_rootPitch7_UMNAST[1,:], tip_u3VsU_rootPitch7_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_u3VsU_rootPitch7_Exp[1,:], tip_u3VsU_rootPitch7_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        scatter!(sqrt.(2*tip_u3Vsq_rootPitch7_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch7_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/PazyWingPitchRange_tipOOP.pdf"))

# Tip twist vs. airspeed for root several pitch angles 
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]", xlims=[0,60], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_twist[i,:], c=colors[i], lw=lw, label=false)
    if θ==5
        plot!(tip_thetaVsU_rootPitch5_UMNAST[1,:], tip_thetaVsU_rootPitch5_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_thetaVsU_rootPitch5_Exp[1,:], tip_thetaVsU_rootPitch5_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
    elseif θ==7
        plot!(tip_thetaVsU_rootPitch7_UMNAST[1,:], tip_thetaVsU_rootPitch7_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_thetaVsU_rootPitch7_Exp[1,:], tip_thetaVsU_rootPitch7_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
    end
end
display(plt2)
savefig(string(absPath,"/PazyWingPitchRange_tipTwist.pdf"))

# Tip AoA vs. airspeed for root several pitch angles 
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Tip AoA [deg]", xlims=[0,60], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_AoA[i,:], c=colors[i], lw=lw, label=false)
end
display(plt3)
savefig(string(absPath,"/PazyWingPitchRange_tipAoA.pdf"))

# Tip in-plane displacement vs. airspeed for root several pitch angles 
plt4 = plot(xlabel="Airspeed [m/s]", ylabel="Tip IP disp. [% semispan]", xlims=[0,60], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_IP[i,:]/L*100, c=colors[i], lw=lw, label=false)
end
display(plt4)
savefig(string(absPath,"/PazyWingPitchRange_tipIP.pdf"))

# Lift coefficient over span
plot_steady_outputs(problem,outputs=["cl","κ2","F3"],colorScheme=:grays,lw=lw,save=true,saveFolder=string(relPath,"/"))

println("Finished PazyWingPitchRangePlotGenerator.jl")