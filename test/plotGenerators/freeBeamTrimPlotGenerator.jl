using Plots

# Run the script
include("../examples/freeBeamTrim.jl")

# Set paths
relPath = "/test/outputs/figures/freeBeamTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,scale=1e2,showScale=true,scalePos=[0.2,0.22],save=true,savePath=string(relPath,"/freeBeamTrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# F3
plt1 = plot(x1/L, F3/F, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3/F\$")
display(plt1)
savefig(string(absPath,"/freeBeamTrim_F.pdf"))

# M2
plt2 = plot(x1/L, M2/(F*L/4), lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2/(FL/4)\$")
display(plt2)
savefig(string(absPath,"/freeBeamTrim_M.pdf"))

println("Finished freeBeamTrimPlotGenerator.jl")