using Plots

# Run the script
include("../examples/midLoadedBeamTrim.jl")

# Set paths
relPath = "/test/outputs/figures/midLoadedBeamTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,scale=1e2,scalePos=[0.1;-0.7;0.05],showScale=true,save=true,savePath=string(relPath,"/midLoadedBeamTrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot(x1/L, u3/(F*L/(48*EI)), lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3 / (FL/48EI)\$")
display(plt1)
savefig(string(absPath,"/midLoadedBeamTrim_u3.pdf"))

# F3
plt2 = plot(x1/L, F3/(F/2), lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3/(F/2)\$ [N]")
display(plt2)
savefig(string(absPath,"/midLoadedBeamTrim_F3.pdf"))

# M2
plt3 = plot(x1/L, M2/(F*L/4), lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2/(FL/4)\$")
display(plt3)
savefig(string(absPath,"/midLoadedBeamTrim_M2.pdf"))

println("Finished midLoadedBeamTrimPlotGenerator.jl")