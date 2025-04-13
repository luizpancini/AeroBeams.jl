using Plots, ColorSchemes

# Run the script
include("../examples/PazyWingFlutterPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingFlutterPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ts = 10
fs = 16
lfs = 9
lw = 2
ms = 10
ms2 = 4
msw = 0
gr()

# Flutter onset and offset speeds vs root pitch angle
humpMode = 3
x1 = [flutterOnsetSpeedsOfMode[i,humpMode][1] for i in eachindex(θRange)]
x2 = [flutterOffsetSpeedsOfMode[i,humpMode][1] for i in eachindex(θRange)]
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]", xlims=[30,90], ylims=[0,7.25], xticks=collect(30:10:90), yticks=collect(0:1:7), legend=:topright, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!(Shape(vcat(x1,reverse(x2)),vcat(θRange,reverse(θRange))), fillcolor = plot_color(:red, 0.25), lw=lw, label="AeroBeams hump flutter region")
scatter!(flutterOnsetVelUp, rootPitchVelUp, shape=:rtriangle, mc=:red, ms=ms, msw=msw, label="Test onset up")
scatter!(flutterOffsetVelUp, rootPitchVelUp, shape=:rtriangle, mc=:green, ms=ms, msw=msw, label="Test offset up")
scatter!(flutterOnsetVelDown, rootPitchVelDown, shape=:ltriangle, mc=:red, ms=ms, msw=msw, label="Test onset down")
scatter!(flutterOffsetVelDown, rootPitchVelDown, shape=:ltriangle, mc=:green, ms=ms, msw=msw, label="Test offset down")
plot!(flutterBoundaryPitch_UMNAST[1,:], flutterBoundaryPitch_UMNAST[2,:], c=:olive, ls=:dash, lw=lw, marker=:circle, ms=ms2, msw=msw, label="UM/NAST (Exp. loss)")
plot!(flutterBoundaryPitch_UMNAST_PanelCoeffs[1,:], flutterBoundaryPitch_UMNAST_PanelCoeffs[2,:], c=:brown, ls=:dash, lw=lw, marker=:circle, ms=ms2, msw=msw, label="UM/NAST (Panel coeffs.)")
plot!(flutterBoundary_UVsPitch_Lambda0_Sharpy[1,1:19], flutterBoundary_UVsPitch_Lambda0_Sharpy[2,1:19], c=:magenta, ls=:dash, lw=lw, marker=:diamond, ms=ms2, msw=msw, label="Sharpy (VLM)")
display(plt1)
savefig(string(absPath,"/PazyWingFlutterPitchRange_flutterBoundaryPitch.pdf"))

# Flutter onset and offset speeds vs tip OOP displacement for varying root pitch angle
θ2plot = [0.5,1,2,3,5,7]
indθ2plot = findall(x -> any(isapprox(x, θ; atol=1e-10) for θ in θ2plot), θRange)
θcolors = cgrad([:blue, :red], length(θ2plot), categorical=true)
xθ = [ 5,  9, 12, 14, 16, 17.5]
yθ = [59, 55, 48, 44, 38, 34]
x1 = [flutterOnsetDispOfMode[i,humpMode][1] for i in eachindex(θRange)]
x2 = [flutterOffsetDispOfMode[i,humpMode][1] for i in eachindex(θRange)]
y1 = [flutterOnsetSpeedsOfMode[i,humpMode][1] for i in eachindex(θRange)]
y2 = [flutterOffsetSpeedsOfMode[i,humpMode][1] for i in eachindex(θRange)]
plt2 = plot(xlabel="Tip OOP displacement [% semispan]", ylabel="Airspeed [m/s]", xlims=[0,32], ylims=[30,90], xticks=collect(0:5:30), yticks=collect(30:10:90), legend=:bottomleft, tickfont=font(ts), guidefont=font(fs), legendfontsize=10)
plot!(Shape(vcat(x1,reverse(x2)),vcat(y1,reverse(y2))), fillcolor = plot_color(:red, 0.25), lw=lw, label=false)
for (n,ind) in enumerate(indθ2plot)
    θ = θ2plot[n]
    plot!(tip_OOP[ind,:]/L*100, URange, c=θcolors[n], ls=:dash, lw=lw, label=false)
    annotate!([xθ[n]],[yθ[n]], text("\$$(round(θ2plot[n],digits=1)) ^\\circ\$", 10, :bottom, θcolors[n]))
end
scatter!(flutterOnsetDispUp, flutterOnsetVelUp, shape=:rtriangle, mc=:red, ms=ms, msw=msw, label=false)
scatter!(flutterOnsetDispDown, flutterOnsetVelDown, shape=:ltriangle, mc=:red, ms=ms, msw=msw, label=false)
plot!(flutterBoundaryDisp_UMNAST[1,:], flutterBoundaryDisp_UMNAST[2,:], c=:olive, ls=:dash, lw=lw, marker=:circle, ms=ms2, msw=msw, label=false)
plot!(flutterBoundaryDisp_UMNAST_PanelCoeffs[1,:], flutterBoundaryDisp_UMNAST_PanelCoeffs[2,:], c=:brown, ls=:dash, lw=lw, marker=:circle, ms=ms2, msw=msw, label=false)
plot!(flutterBoundary_dispVsU_Lambda0_Sharpy[1,:]*100,flutterBoundary_dispVsU_Lambda0_Sharpy[2,:], c=:magenta, ls=:dash, lw=lw, marker=:diamond, ms=ms2, msw=msw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingFlutterPitchRange_flutterBoundaryDisp.pdf"))

println("Finished PazyWingFlutterPitchRangePlotGenerator.jl")