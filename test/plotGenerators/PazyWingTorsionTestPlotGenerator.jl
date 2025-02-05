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
ts = 10
fs = 16
lw = 2
ms = 4
msw = 0
gr()

# Tip midchord OOP displacement (offset from zero tip mass value) vs. tip mass
plt1 = plot(xlabel="Tip mass [kg]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,3], tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
plot!(mRange, (tip_OOP.-tip_OOP[1])/L*100, c=:black, lw=lw, label="AeroBeams")
plot!(u3UMNAST[1,:], u3UMNAST[2,:], c=:black, ls=:dot, lw=lw, label="UM/NAST - Riso & Cesnik (2022)")
scatter!(u3Exp[1,:], u3Exp[2,:], mc=:black, ms=ms, msw=msw, label="Exp. - Avin et al. (2022)")
display(plt1)
savefig(string(absPath,"/PazyWingTorsionTest_OOP.pdf"))

# Tip twist vs. tip mass
plt2 = plot(xlabel="Tip mass [kg]", ylabel="Tip twist [deg]", xlims=[0,3], tickfont=font(ts), guidefont=font(fs))
plot!(mRange, tip_twist, c=:black, lw=lw, label=false)
plot!(θUMNAST[1,:], θUMNAST[2,:], c=:black, ls=:dot, lw=lw, label=false)
scatter!(θExp[1,:], θExp[2,:], mc=:black, ms=ms, msw=msw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingTorsionTest_twist.pdf"))

println("Finished PazyWingTorsionTestPlotGenerator.jl")