using Plots

# Run the script
include("../examples/PazyWingContinuous2DSpaceGust.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingContinuous2DSpaceGust"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=10,plotLimits=[(-L/2,L/2),(-L/2,L/2),(0,L)],save=true,savePath=string(relPath,"/PazyWingContinuous2DSpaceGust_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
gr()

# Tip displacement
plt1 = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]")
plot!(t, tipOOP/L*100, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/PazyWingContinuous2DSpaceGust_disp.pdf"))

# Tip AoA
plt2 = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]")
plot!(t, tipAoA*180/π, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingContinuous2DSpaceGust_AoA.pdf"))

# 3/4-span cn
plt3 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_n\$")
plot!(t, tqSpan_cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/PazyWingContinuous2DSpaceGust_cn.pdf"))

# 3/4-span cm
plt4 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_m\$")
plot!(t, tqSpan_cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/PazyWingContinuous2DSpaceGust_cm.pdf"))

# 3/4-span ct
plt5 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_t\$")
plot!(t, tqSpan_ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/PazyWingContinuous2DSpaceGust_ct.pdf"))

# Aero states at 3/4-span
nAeroStates = problem.model.elements[1].aero.nTotalAeroStates
colors = get(colorschemes[:rainbow], LinRange(0, 1, nAeroStates))
tqsχ_ = Array{Vector{Float64}}(undef,nAeroStates)
for i in 1:nAeroStates
    tqsχ_[i] = [tqsχ[tt][i] for tt in 1:length(t)]
end
plt6 = plot(xlabel="Time [s]", ylabel="")
for i in 1:nAeroStates
    plot!(t, tqsχ_[i], c=colors[i], lw=lw, label="\$\\chi $(i)\$")
end
display(plt6)
savefig(string(absPath,"/PazyWingContinuous2DSpaceGust_states.pdf"))

# "Vertical" gust velocity
W = t -> t0<t<t0+τ ? gust.W([0; t*U; 0]) : 0
plt7 = plot(xlabel="Time [s]", ylabel="Gust velocity [m/s]")
plot!(t, W.(t), color=:black, lw=lw, label=false)
display(plt7)
savefig(string(absPath,"/PazyWingContinuous2DSpaceGust_gust.pdf"))

println("Finished PazyWingContinuous2DSpaceGustPlotGenerator.jl")