using Plots, ColorSchemes

# Run the script
include("../examples/elbowFrame.jl")

# Set paths
relPath = "/test/outputs/figures/elbowFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,plotFrequency=5,view=(30,30),plotLimits=[(0,L),(0,L),(-L/2,L/2)],save=true,savePath=string(relPath,"/elbowFrame_deformation.gif"),displayProgress=true)

# Plot configurations
y = [u3_elbow, u3_tip]
labels = ["Elbow" "Tip"]
colors = [:blue,:orange]
lw = 2
ms = 4
msw = 0
gr()

# Plot displacements over time
plt1 = plot( xlabel="\$t\$ [s]", ylabel="\$u_3\$ [in]", title="OOP displacements",legend=:bottomleft)
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Simo & Vu-Quoc (1987)")
for i=1:2
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(t, y, lw=lw, palette=colors, label=false)
scatter!([u3ElbowRef[1,:],u3TipRef[1,:]], [u3ElbowRef[2,:],u3TipRef[2,:]], palette=colors,ms=ms,msw=msw,label=false)
display(plt1)
savefig(string(absPath,"/elbowFrame_disp.pdf"))

println("Finished elbowFrame.jl")