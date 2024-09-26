using Plots, ColorSchemes

# Run the script
include("../examples/SMWFlutterPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/SMWFlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
flutterModes = [2,5]
flutterOnsetModesLabels = ["2nd OOP - onset" "2nd T-IP - onset"]
flutterOffsetModesLabels = ["2nd OOP - offset" "2nd T-IP - offset"]
gr()

# Flutter speed and frequency vs. root pitch angle
plt11 = plot(ylabel="Flutter speed [m/s]", xlims=[0,5], ylims=[0,35], legend=:bottomleft)
for (i,m) in enumerate(flutterModes)
    plot!(θRange, flutterOnsetSpeed[:,m], c=modeColors[m], ls=:solid, lw=lw, label=flutterOnsetModesLabels[i])
    plot!(θRange, flutterOffsetSpeed[:,m], c=modeColors[m], ls=:dash, lw=lw, label=flutterOffsetModesLabels[i])
end
plot!(flutterSpeedRef[1,:], flutterSpeedRef[2,:], c=:black, ls=:solid, lw=lw, label="Patil et al. (2001)")
plt12 = plot(xlabel="Root pitch angle, \$\\theta\$ [deg]", ylabel="Flutter frequency [rad/s]", xlims=[0,5], ylims=[0,35])
for (i,m) in enumerate(flutterModes)
    plot!(θRange, flutterOnsetFreq[:,m], c=modeColors[m], ls=:solid, lw=lw, label=false)
    plot!(θRange, flutterOffsetFreq[:,m], c=modeColors[m], ls=:dash, lw=lw, label=false)
end
plot!(flutterFreqRef[1,:], flutterFreqRef[2,:], c=:black, ls=:solid, lw=lw, label=false)
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/SMWFlutterPitchRange_flutterSpeed.pdf"))

# Flutter tip displacement vs. root pitch angle
plt2 = plot(xlabel="Root pitch angle, \$\\theta\$ [deg]",ylabel="Flutter tip displacement [m]", xlims=[0,5], ylims=[-3,3])
for (i,m) in enumerate(flutterModes)
    plot!(θRange, flutterOnsetTipDisp[:,m], c=modeColors[m], ls=:solid, lw=lw, label=false)
    plot!(θRange, flutterOffsetTipDisp[:,m], c=modeColors[m], ls=:dash, lw=lw, label=false)
end
plot!(flutterTipDispRef[1,:], flutterTipDispRef[2,:], c=:black, ls=:solid, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/SMWFlutterPitchRange_flutterDisp.pdf"))

# Tip displacement vs airspeed
θplot = [0; 0.5; 1.0; 2.0]
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(θplot)))
indθ2plot = findall(vec(any(θRange .== θplot', dims=2)))
plt3 = plot(xlabel="Tip displacement [m]",ylabel="Airspeed [m/s]", xlims=[-3,3], ylims=[0,35])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="Patil et al. (2001)")
for i in eachindex(indθ2plot)
    plot!(tip_u3[indθ2plot[i],:], URange, c=colors[i], lw=lw, label="\\theta=$(θplot[i]) deg")
    if θplot[i] == 0.0
        scatter!(speedVsDispRootAoA0[1,:],speedVsDispRootAoA0[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 0.5
        scatter!(speedVsDispRootAoA05[1,:],speedVsDispRootAoA05[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 1.0
        scatter!(speedVsDispRootAoA1[1,:],speedVsDispRootAoA1[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 2.0
        scatter!(speedVsDispRootAoA2[1,:],speedVsDispRootAoA2[2,:], mc=colors[i], ms=ms, msw=0, label=false)    
    end
end
display(plt3)
savefig(string(absPath,"/SMWFlutterPitchRange_tipDisp.pdf"))

# V-g-f of selected mode for selected root angles
mode2plot = 2
θplot = [2.0; 3.0; 5.0]
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(θplot)))
indθ2plot = findall(vec(any(θRange .== θplot', dims=2)))
plt41 = plot(ylabel="Frequency [rad/s]", xlims=[10,30], legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="Patil et al. (2001)")
for i in eachindex(indθ2plot)
    plot!(URange, [freqs[indθ2plot[i],j][mode2plot] for j in eachindex(URange)], c=colors[i], lw=lw, label=false)
    if θplot[i] == 2.0
        scatter!(freqVsSpeedRootAoA2[1,:],freqVsSpeedRootAoA2[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 3.0
        scatter!(freqVsSpeedRootAoA3[1,:],freqVsSpeedRootAoA3[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 5.0
        scatter!(freqVsSpeedRootAoA5[1,:],freqVsSpeedRootAoA5[2,:], mc=colors[i], ms=ms, msw=0, label=false)   
    end
end
plt42 = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", xlims=[10,30], ylims=[-1.,0.25], legend=:bottomleft)
for i in eachindex(indθ2plot)
    plot!(URange, [damps[indθ2plot[i],j][mode2plot] for j in eachindex(URange)], c=colors[i], lw=lw, label="\\theta=$(θplot[i]) deg")
    if θplot[i] == 2.0
        scatter!(dampVsSpeedRootAoA2[1,:],dampVsSpeedRootAoA2[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 3.0
        scatter!(dampVsSpeedRootAoA3[1,:],dampVsSpeedRootAoA3[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 5.0
        scatter!(dampVsSpeedRootAoA5[1,:],dampVsSpeedRootAoA5[2,:], mc=colors[i], ms=ms, msw=0, label=false)   
    end
end
plt4 = plot(plt41,plt42, layout=(2,1))
display(plt4)
savefig(string(absPath,"/SMWFlutterPitchRange_Vgf.pdf"))

println("Finished SMWFlutterPitchRangePlotGenerator.jl")