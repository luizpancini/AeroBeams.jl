using Plots

# Run the script
include("../examples/HealyFFWTOMCGustFloating.jl")

# Set paths
relPath = "/test/outputs/figures/HealyFFWTOMCGustFloating"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Beam geometry
L = HealyFFWTOMCGustFloating.beams[1].length

# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=20,fps=30,view=(30,30),plotLimits=([0,L],[-L/2,L/2],[-L/2,L/2]),plotDistLoads=false,save=true,savePath=string(relPath,"/HealyFFWTOMCGustFloating_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
gr()

# Steady p1
plt_p1 = plot(xlabel="\$x_1/L\$", ylabel="Steady \$p_1\$")
plot!(x1/L, p1, c=:black, lw=lw, label=false)
display(plt_p1)
savefig(string(absPath,"/HealyFFWTOMCGustFloating_p1.pdf"))

# Steady p2
plt_p2 = plot(xlabel="\$x_1/L\$", ylabel="Steady \$p_2\$")
plot!(x1/L, p2, c=:black, lw=lw, label=false)
display(plt_p2)
savefig(string(absPath,"/HealyFFWTOMCGustFloating_p2.pdf"))

# Steady p3
plt_p3 = plot(xlabel="\$x_1/L\$", ylabel="Steady \$p_3\$")
plot!(x1/L, p3, c=:black, lw=lw, label=false)
display(plt_p3)
savefig(string(absPath,"/HealyFFWTOMCGustFloating_p3.pdf"))

# Steady cn over span
plt_scn = plot(xlabel="\$x_1/L\$", ylabel="Steady \$c_n\$")
plot!(x1_e/L, steady_cn_over_span, c=:black, lw=lw, label=false)
display(plt_scn)
savefig(string(absPath,"/HealyFFWTOMCGustFloating_scn.pdf"))

# Tip displacement
plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [m]")
plot!(t, tipOOP, c=:black, lw=lw, label=false)
display(plt_tipOOP)
savefig(string(absPath,"/HealyFFWTOMCGustFloating_tipdisp.pdf"))

# Root OOP bending moment increment
plt_ΔM2root = plot(xlabel="Time [s]", ylabel="ΔWRBM [N.m]")
plot!(t, -(M2_root .- M2_root[1]), c=:black, lw=lw, label=false)
display(plt_ΔM2root)
savefig(string(absPath,"/HealyFFWTOMCGustFloating_DeltaM2root.pdf"))

# AoA
plt_aoa = plot(xlabel="Time [s]", ylabel="Angle of attack [deg]")
plot!(t, root_αₑ*180/π, lw=lw, label="Root")
plot!(t, tqSpan_αₑ*180/π, lw=lw, label="3/4-span")
plot!(t, tip_αₑ*180/π, lw=lw, label="Tip")
display(plt_aoa)
savefig(string(absPath,"/HealyFFWTOMCGustFloating_AoA.pdf"))

# cn
plt_cn = plot(xlabel="Time [s]", ylabel="\$c_n\$")
plot!(t, root_cn, lw=lw, label="Root")
plot!(t, tqSpan_cn, lw=lw, label="3/4-span")
plot!(t, tip_cn, lw=lw, label="Tip")
display(plt_cn)
savefig(string(absPath,"/HealyFFWTOMCGustFloating_cn.pdf"))

# cm
plt_cm = plot(xlabel="Time [s]", ylabel="\$c_m\$")
plot!(t, root_cm, lw=lw, label="Root")
plot!(t, tqSpan_cm, lw=lw, label="3/4-span")
plot!(t, tip_cm, lw=lw, label="Tip")
display(plt_cm)
savefig(string(absPath,"/HealyFFWTOMCGustFloating_cm.pdf"))

# 3/4-span ct
plt_ct = plot(xlabel="Time [s]", ylabel="\$c_t\$")
plot!(t, root_ct, lw=lw, label="Root")
plot!(t, tqSpan_ct, lw=lw, label="3/4-span")
plot!(t, tip_ct, lw=lw, label="Tip")
display(plt_ct)
savefig(string(absPath,"/HealyFFWTOMCGustFloating_ct.pdf"))

# Aero states at 3/4-span
nAeroStates = problem.model.elements[1].aero.nTotalAeroStates
colors = get(colorschemes[:rainbow], LinRange(0, 1, nAeroStates))
tqsχ_ = Array{Vector{Float64}}(undef,nAeroStates)
for i in 1:nAeroStates
    tqsχ_[i] = [tqsχ[tt][i] for tt in 1:length(t)]
end
plt_aeroStates = plot(xlabel="Time [s]", ylabel="")
for i in 1:nAeroStates
    plot!(t, tqsχ_[i], c=colors[i], lw=lw, label="\$\\chi $(i)\$")
end
display(plt_aeroStates)
savefig(string(absPath,"/HealyFFWTOMCGustFloating_states.pdf"))

println("Finished HealyFFWTOMCGustFloatingPlotGenerator.jl")