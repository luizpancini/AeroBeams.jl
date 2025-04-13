using Plots

# Run the script
include("../examples/plungingAirfoilDSModelTest.jl")

# Set paths
relPath = "/test/outputs/figures/plungingAirfoilDSModelTest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ts = 10
fs = 16
lw = 2
ms = 3
msw = 0
gr()

# Index range of last cycle
rangeLastCycle = ceil(Int,(tf-T)/Δt):Nt

# cn vs α_eq
plt_cn = plot(xlabel="\$\\alpha_{eq}\$ [deg]", ylabel="\$c_n\$", tickfont=font(ts), guidefont=font(fs))
plot!(α_eq[rangeLastCycle]*180/π, cn[rangeLastCycle], c=:black, ls=:solid, lw=lw, label=false)
plot!(cn_num[1,:], cn_num[2,:], c=:black, ls=:dashdot, lw=lw, label=false)
scatter!(cn_exp[1,:], cn_exp[2,:], c=:black, ms=ms, msw=msw, label=false)
display(plt_cn)
savefig(string(absPath,"/plungingAirfoilDSModelTest_cn.pdf"))

# cm vs α_eq
plt_cm = plot(xlabel="\$\\alpha_{eq}\$ [deg]", ylabel="\$c_m\$", tickfont=font(ts),  guidefont=font(fs), legend=:bottomleft, legendfontsize=12)
plot!(α_eq[rangeLastCycle]*180/π, cm[rangeLastCycle], c=:black, ls=:solid, lw=lw, label="AeroBeams")
plot!(cm_num[1,:], cm_num[2,:], c=:black, ls=:dashdot, lw=lw, label="Model - Tyler & Leishman (1992)")
scatter!(cm_exp[1,:], cm_exp[2,:], c=:black, ms=ms, msw=msw, label="Exp. - Liiva et al. (1968)")
display(plt_cm)
savefig(string(absPath,"/plungingAirfoilDSModelTest_cm.pdf"))

println("Finished plungingAirfoilDSModelTestPlotGenerator.jl")