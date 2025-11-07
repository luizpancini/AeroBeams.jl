using AeroBeams, Plots

# Indicial function
circulatoryIndicialFunction = "Wagner"
solver = Indicial(circulatoryIndicialFunction=circulatoryIndicialFunction)
τ = 0:1e-2:20
Φ = τ -> 1 - sum([solver.AC[i]*exp(-solver.bC[i]*τ) for i in 1:solver.nCirculatoryStates])
Φnum = Φ.(τ)

# Set paths
relPath = "/dev/misc/outputs/figures/indicialFunction"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ts = 20
fs = 30
lw = 4
fps = 30
pltFreq = 20
DPI = 300

# Animated plot of curvature
gr()
τSpaced = τ[1:pltFreq:end]
ΦnumSpaced = Φnum[1:pltFreq:end]
plt = plot(xlabel="\$\\tau\$", ylabel="\$\\Phi\$", title="\$\\Phi=\\Delta\\alpha^E/\\Delta\\alpha^{QS}\$", xlims=[0,τ[end]], ylims=[0, 1], xticks=false, tickfont=font(ts), guidefont=font(fs), titlefont=font(fs), dpi=DPI)
anim = @animate for (tt,timeNow) in enumerate(τSpaced)
    plot!(τSpaced[1:tt], ΦnumSpaced[1:tt], c=:forestgreen, lw=lw, label=false)
end
gif_handle = gif(anim, string(absPath,"/indicialFunction.gif"), fps=fps)
display(gif_handle)