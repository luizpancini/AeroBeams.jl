using Plots

# Run the script
# include("../examples/PazyWingFlutterPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingFlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 10
msw = 0
gr()

# Flutter onset and offset speeds vs root pitch angle
mode2plot = 3
x1 = [flutterOnsetSpeedsOfMode[i,mode2plot][1] for i in eachindex(θRange)]
x2 = [flutterOffsetSpeedsOfMode[i,mode2plot][1] for i in eachindex(θRange)]
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]", xlims=[0,90], ylims=[0,7.25], xticks=collect(0:15:90), yticks=collect(0:1:7), legend=:bottomleft)
plot!(Shape(vcat(x1,reverse(x2[3:end]),97,100),vcat(θRange,reverse(θRange))), fillcolor = plot_color(:red, 0.25), lw=lw, label="AeroBeams flutter region")
scatter!(flutterOnsetVelUp, rootPitchVelUp, shape=:rtriangle, mc=:red, ms=ms, msw=msw, label="Test onset up")
scatter!(flutterOffsetVelUp, rootPitchVelUp, shape=:rtriangle, mc=:green, ms=ms, msw=msw, label="Test offset up")
scatter!(flutterOnsetVelDown, rootPitchVelDown, shape=:ltriangle, mc=:red, ms=ms, msw=msw, label="Test onset down")
scatter!(flutterOffsetVelDown, rootPitchVelDown, shape=:ltriangle, mc=:green, ms=ms, msw=msw, label="Test offset down")
display(plt1)
savefig(string(absPath,"/PazyWingFlutterPitchRange_flutterBoundaryPitch.pdf"))

# Flutter onset and offset speeds vs tip OOP displacement for varying root pitch angle
θ2plot = [0.5,1,2,3,5,7]
indθ2plot = findall(vec(any(θRange .== θ2plot', dims=2)))
x1 = [flutterOnsetDispOfMode[i,mode2plot][1] for i in eachindex(θRange)]
x2 = [flutterOffsetDispOfMode[i,mode2plot][1] for i in eachindex(θRange)]
y1 = [flutterOnsetSpeedsOfMode[i,mode2plot][1] for i in eachindex(θRange)]
y2 = [flutterOffsetSpeedsOfMode[i,mode2plot][1] for i in eachindex(θRange)]
plt2 = plot(xlabel="Tip OOP displacement [% semispan]", ylabel="Airspeed [m/s]", xlims=[0,32], ylims=[30,90], xticks=collect(0:5:30), yticks=collect(30:10:90), legend=:bottomleft)
plot!(Shape(vcat(x1,reverse(x2[3:end]),22,22,0),vcat(y1,reverse(y2[3:end]),97,100,90)), fillcolor = plot_color(:red, 0.25), lw=lw, label="AeroBeams flutter region")
for (n,ind) in enumerate(indθ2plot)
    θ = θ2plot[n]
    plot!(tip_OOP[ind,:]/L*100, URange, c=:black, ls=:dash, lw=lw, label=false)
    if θ==0.5
        xind,yind = 16,79
    elseif θ==1
        xind,yind = 20,72
    elseif θ==2
        xind,yind = 22.5,63
    elseif θ==3
        xind,yind = 24,56
    elseif θ==5
        xind,yind = 26,49
    elseif θ==7
        xind,yind = 27,42.5  
    end
    annotate!([xind],[yind], text("$(θ2plot[n]) deg", 10, :bottom))
end
scatter!(flutterOnsetDispUp, flutterOnsetVelUp, shape=:rtriangle, mc=:red, ms=ms, msw=msw, label="Test onset up")
scatter!(flutterOnsetDispDown, flutterOnsetVelDown, shape=:ltriangle, mc=:red, ms=ms, msw=msw, label="Test onset down")
display(plt2)
savefig(string(absPath,"/PazyWingFlutterPitchRange_flutterBoundaryDisp.pdf"))

println("Finished PazyWingFlutterPitchRangePlotGenerator.jl")