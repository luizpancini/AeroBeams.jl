using Plots, ColorSchemes

# Run the script
include("../examples/SMWFlutterPrecurvatureRange2.jl")

# Set paths
relPath = "/test/outputs/figures/SMWFlutterPrecurvatureRange2"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(kRange)))
lw = 2
ms = 3
gr()

# Initialize arrays
fSpeedOnsetMode = Array{Vector{Float64}}(undef,nModes,length(kRange))
fSpeedOffsetMode = Array{Vector{Float64}}(undef,nModes,length(kRange))

# Flutter speed vs. tip load
plt1 = plot(xlabel="Root angle [deg]", ylabel="Flutter speed [m/s]", xlims=[0,5], ylims=[0,40])
for (i,k) in enumerate(kRange)
    plot!([NaN], [NaN], c=colors[i], ls=:solid, lw=lw, label="k=$k - onset")
    plot!([NaN], [NaN], c=colors[i], ls=:dash, lw=lw, label="k=$k - offset")
end
for m in 1:nModes
    for (i,k) in enumerate(kRange)
        fSpeedOnsetMode[m,i] = vcat([flutterOnsetSpeed[i][j][m][1] for j in eachindex(θRange)]...)
        fSpeedOffsetMode[m,i] = vcat([flutterOffsetSpeed[i][j][m][1] for j in eachindex(θRange)]...)
        plot!(θRange, fSpeedOnsetMode[m,i], c=colors[i], ls=:solid, lw=lw, label=false)
        plot!(θRange, fSpeedOffsetMode[m,i], c=colors[i], ls=:dash, lw=lw, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/SMWFlutterPrecurvatureRange2_flutter.pdf"))

println("Finished SMWFlutterPrecurvatureRange2PlotGenerator.jl")