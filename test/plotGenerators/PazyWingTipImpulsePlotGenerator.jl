using Plots, ColorSchemes

# Run the script
include("../examples/PazyWingTipImpulse.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingTipImpulse"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,refBasis="A",plotFrequency=50,plotLimits=([-L/2,L/2],[-L/2,L/2],[0,L]),save=true,savePath=string(relPath,"/PazyWingTipImpulse_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Tip displacement
plt1 = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]")
plot!(t, tipOOP/L*100, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/PazyWingTipImpulse_disp.pdf"))

# Tip AoA
plt2 = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]")
plot!(t, tipAoA*180/π, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingTipImpulse_AoA.pdf"))

# 3/4-span cn
plt3 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_n\$")
plot!(t, tqSpan_cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/PazyWingTipImpulse_cn.pdf"))

# 3/4-span cm
plt4 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_m\$")
plot!(t, tqSpan_cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/PazyWingTipImpulse_cm.pdf"))

# 3/4-span ct
plt5 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_t\$")
plot!(t, tqSpan_ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/PazyWingTipImpulse_ct.pdf"))

# Aero states at 3/4-span
nAeroStates = problem.model.elements[1].aero.nTotalAeroStates
colors = get(colorschemes[:rainbow], LinRange(0, 1, 8))
tqsχ_ = Array{Vector{Float64}}(undef,nAeroStates)
for i in 1:nAeroStates
    tqsχ_[i] = [tqsχ[tt][i] for tt in 1:length(t)]
end
plt6 = plot(xlabel="Time [s]", ylabel="")
for i in 1:nAeroStates
    plot!(t, tqsχ_[i], c=colors[i], lw=lw, label="\$\\chi $(i)\$")
end
display(plt6)
savefig(string(absPath,"/PazyWingTipImpulse_states.pdf"))

# Aero states' rates at 3/4-span
tqsχdot_ = Array{Vector{Float64}}(undef,nAeroStates)
for i in 1:nAeroStates
    tqsχdot_[i] = [tqsχdot[tt][i] for tt in 1:length(t)]
end
plt7 = plot(xlabel="Time [s]", ylabel="")
for i in 1:nAeroStates
    plot!(t, tqsχdot_[i], c=colors[i], lw=lw, label="\$\\dot{\\chi} $(i)\$")
end
display(plt7)
savefig(string(absPath,"/PazyWingTipImpulse_statesRates.pdf"))

# Root, 3/4-span and tip lift and drag coefficients over time
plot_time_outputs(problem,elements=[1,12,nElem],elementalOutputs=["cl","cd"],save=true,saveFolder=string(relPath,"/"))

println("Finished PazyWingTipImpulse.jl")