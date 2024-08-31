using Plots

# Run the script
include("../examples/PazyWingTorsionTest.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingTorsionTest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,view=(30,30),save=true,savePath=string(relPath,"/PazyWingTorsionTest_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
ms = 4
msw = 0
gr()

# Tip midchord OOP displacement (offset from zero tip mass value) vs. tip mass
plt1 = plot(xlabel="Tip mass [kg]", ylabel="Tip OOP displacement offset [% semispan]", xlims=[0,3])
plot!(mRange, (tip_OOP.-tip_OOP[1])/L*100, c=:black, lw=lw, label="AeroBeams")
plot!(u3UMNAST[1,:], u3UMNAST[2,:], c=:blue, ls=:dash, lw=lw, label="UM/NAST")
scatter!(u3Exp[1,:], u3Exp[2,:], mc=:red, ms=ms, msw=msw, label="Exp.")
display(plt1)
savefig(string(absPath,"/PazyWingTorsionTest_OOP.pdf"))

# Tip twist vs. tip mass
plt2 = plot(xlabel="Tip mass [kg]", ylabel="Tip twist [deg]", xlims=[0,3])
plot!(mRange, tip_twist, c=:black, lw=lw, label="AeroBeams")
plot!(θUMNAST[1,:], θUMNAST[2,:], c=:blue, ls=:dash, lw=lw, label="UM/NAST")
scatter!(θExp[1,:], θExp[2,:], mc=:red, ms=ms, msw=msw, label="Exp.")
display(plt2)
savefig(string(absPath,"/PazyWingTorsionTest_twist.pdf"))

println("Finished PazyWingTorsionTestPlotGenerator.jl")