using Plots

# Run the script
include("../examples/sweptTipRotor.jl")

# Set paths
relPath = "/test/outputs/figures/sweptTipRotor"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes
modesPlot = plot_mode_shapes(problem,scale=5,view=(30,30),legendPos=:best,frequencyLabel="frequency",save=true,savePath=string(relPath,"/sweptTipRotor_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
gr()
colors = get(colorschemes[:darkrainbow], LinRange(0, 1, length(ωRange)))
lgdtitle = "Lines: Numerical\nMarkers: Exp. - Epps & Chandra (1996)"
lw = 2
ms = 5
msw = 0

# Plot 1st bending mode frequency over tip angle for several angular velocities
plt1 = plot(xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]")
for (i,ω) in enumerate(ωRange) 
    mode = 1
    ωRPM = round(Int,ω/(2*π/60)) 
    numFreqs1 = [numFreqs[i,j][mode]/(2π) for j in 1:size(numFreqs, 2)]
    plot!(tipAngleRange*180/π,numFreqs1, lc=colors[i], lw=lw, label=false)
    scatter!(expTipAngles,expFreqs1[i,:], mc=colors[i], ms=ms, msw=msw, label=false)
    plot!([NaN], [NaN], lc=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label="\$\\omega\$ = $ωRPM rpm")
end
plot!(title="1st bending")
# plot!(legendtitle=lgdtitle)
display(plt1)
savefig(string(absPath,"/sweptTipRotor_1B.pdf"))

# Plot 2nd bending mode frequency over tip angle for several angular velocities
plt2 = plot(xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]")
for (i,ω) in enumerate(ωRange) 
    if i < 3
        mode = 2
    else
        mode = 3
    end
    ωRPM = round(Int,ω/(2*π/60)) 
    numFreqs2 = [numFreqs[i,j][mode]/(2π) for j in 1:size(numFreqs, 2)]
    plot!(tipAngleRange*180/π,numFreqs2, lc=colors[i], lw=lw, label=false)
    scatter!(expTipAngles,expFreqs2[i,:], mc=colors[i], ms=ms, msw=msw, label=false)
    plot!([NaN], [NaN], lc=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label="\$\\omega\$ = $ωRPM rpm")
end
plot!(title="2nd bending")
# plot!(legendtitle=lgdtitle)
display(plt2)
savefig(string(absPath,"/sweptTipRotor_2B.pdf"))

# Plot 3rd bending mode frequency over tip angle for several angular velocities
plt3 = plot(xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]")
for (i,ω) in enumerate(ωRange) 
    mode = 4
    ωRPM = round(Int,ω/(2*π/60)) 
    numFreqs3 = [numFreqs[i,j][mode]/(2π) for j in 1:size(numFreqs, 2)]
    plot!(tipAngleRange*180/π,numFreqs3, lc=colors[i], lw=lw, label=false)
    scatter!(expTipAngles,expFreqs3[i,:], mc=colors[i], ms=ms, msw=msw, label=false)
    plot!([NaN], [NaN], lc=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label="\$\\omega\$ = $ωRPM rpm")
end
plot!(title="3rd bending")
# plot!(legendtitle=lgdtitle)
display(plt3)
savefig(string(absPath,"/sweptTipRotor_3B.pdf"))

# Plot coupled bending-torsion modes frequency over tip angle for ω = 750 rpm
plt4 = plot()
modes = [5,6,7]
modeLabels = ["1T/5B","5B/1T","4B/1T"]
for (i,ω) in enumerate(ωRange) 
    ωRPM = round(Int,ω/(2π/60)) 
    if i == 3
        numFreqsCoupled = zeros(length(tipAngleRange))
        for j in 1:size(numFreqs, 2)
            if tipAngleRange[j]*180/π >= 20 
                mode = 8
            else
                mode = 7
            end
            numFreqsCoupled[j] = numFreqs[end,j][mode]/(2*π)
        end
    else
        mode = modes[i]
        numFreqsCoupled = [numFreqs[end,j][mode]/(2*π) for j in 1:size(numFreqs, 2)]
    end
    plot!(tipAngleRange*180/π,numFreqsCoupled, lc=colors[i], lw=lw, xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]", label=false)
    scatter!(expTipAngles,expFreqs4[i,:], mc=colors[i], ms=ms, msw=msw, label=false)
    plot!([NaN], [NaN], lc=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=modeLabels[i])
end
plot!(title="Coupled bending-torsion at \$\\omega\$ = 750 rpm")
# plot!(legendtitle=lgdtitle)
display(plt4)
savefig(string(absPath,"/sweptTipRotor_TB.pdf"))

println("Finished sweptTipRotorPlotGenerator.jl")