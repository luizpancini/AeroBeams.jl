using Plots, ColorSchemes

# Run the script
include("../examples/SMWFlutterPrecurvatureRange.jl")

# Set paths
relPath = "/test/outputs/figures/SMWFlutterPrecurvatureRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(kRange)))
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
mshape = [:utriangle, :utriangle, :circle, :circle, :utriangle]
lw = 2
ms = 4
gr()

# Flutter speed vs. tip load
plt1 = plot(xlabel="Tip load [N]", ylabel="Flutter speed [m/s]", xlims=[0,40], ylims=[0,40])
for (i,k) in enumerate(kRange)
    plot!(F3Range, flutterSpeed[i,:], c=colors[i], lw=lw, label="AeroBeams - \$k_2\$=$k")
end
plot!(flutterSpeedVsTipLoadk0[1,:], flutterSpeedVsTipLoadk0[2,:], c=:black, ls=:dash, lw=lw, label="Patil et al. (2001)")
plot!(flutterSpeedVsTipLoadk2[1,:], flutterSpeedVsTipLoadk2[2,:], c=:black, ls=:dash, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/SMWFlutterPrecurvatureRange_flutter.pdf"))

# Flutter speed vs. tip disp
plt2 = plot(xlabel="Tip OOP position [% semispan]", ylabel="Flutter speed [m/s]", xlims=[-20,20], ylims=[0,40])
for (i,k) in enumerate(kRange)
    plot!((flutterTipDisp[i,:] .+ SMWFlutterPrecurvatureRange[i,1,1].elements[end].r_n2[3])/16*100, flutterSpeed[i,:], c=colors[i], lw=lw, label="\$k_2\$=$k")
end
display(plt2)
savefig(string(absPath,"/SMWFlutterPrecurvatureRange_flutterTipDisp.pdf"))

# Frequency vs tip displacement of initially straight wing (at lowest airspeed)
plt3 = plot(xlabel = "Tip OOP displacement [m]", ylabel="Frequency [rad/s]", ylims=[0,50])
scatter!([NaN], [NaN], c=:black, shape=mshape[1], ms=ms, msw=msw, label="OOP bending modes")
scatter!([NaN], [NaN], c=:black, shape=mshape[3], ms=ms, msw=msw, label="T-IP bending modes")
for (i,F3) in enumerate(F3Range)
    for mode in 1:nModes
        # Fix mode swap
        if mode == 4 && i > 8
            colorNow = modeColors[5]
            shapeNow = mshape[5]
        elseif mode == 5 && i > 8
            colorNow = modeColors[4]
            shapeNow = mshape[4]
        else
            colorNow=modeColors[mode] 
            shapeNow = mshape[mode]           
        end
        # Plot
        scatter!(tip_u3[1,i,:], [modeFrequencies[1,i,mode][1]], c=colorNow, shape=shapeNow, ms=ms, msw=msw, label=false)
    end
end
display(plt3)
savefig(string(absPath,"/SMWFlutterPrecurvatureRange_freq_vs_disp.pdf"))

println("Finished SMWFlutterPrecurvatureRangePlotGenerator.jl")