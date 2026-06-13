using Plots, ColorSchemes

# Run the script
include("../examples/sweptBeamStaticLoading.jl")

# Set paths
relPath = "/test/outputs/figures/sweptBeamStaticLoading"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 4
gr()

# Plot displaced beams
plt_beams = plot_steady_deformations(problem,backendSymbol=:plotlyjs,plotBCs=true,plotUndeformed=false,plotDistLoads=false,plotAxes=true,plotGrid=true,legendEntries=["\$Λ = $(round(Int,Λ*180/π)) ^\\circ \$" for Λ in ΛRange],legendPos=:right,view=(80,30),plotLimits=([0,L],[-L*3/4,L/4],[-L/2,L/2]),save=true,savePath=string(relPath,"/beam_shapes.pdf"))
display(plt_beams)

# Tip displacement vs. sweep angle
gr()
plt_u3 = plot(xlabel="Sweep angle [deg]", ylabel="Tip \$u_3/L\$")
for i in eachindex(ΛRange)
    plot!(ΛRange*180/π, tip_u3/L, c=:black, lw=lw, label=false)
end
display(plt_u3)
savefig(string(absPath,"/sweptBeamStaticLoadingTipDisp.pdf"))

# Tip twist vs. sweep angle
plt_twist = plot(xlabel="Sweep angle [deg]", ylabel="Tip twist [deg]")
for i in eachindex(ΛRange)
    plot!(ΛRange*180/π, tip_twist, c=:black, lw=lw, label=false)
end
display(plt_twist)
savefig(string(absPath,"/sweptBeamStaticLoadingTipTwist.pdf"))

println("Finished sweptBeamStaticLoadingPlotGenerator.jl")