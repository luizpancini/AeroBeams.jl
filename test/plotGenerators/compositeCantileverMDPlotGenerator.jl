using Plots

# Run the script
include("../examples/compositeCantileverMD.jl")

# Set paths
relPath = "/test/outputs/figures/compositeCantileverMD"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
labels = ["\$-u_1\$" "\$u_2\$" "\$-u_3\$"]
xLabel = "\$-u_1, u_2, -u_3\$ [m]"
yLabel = "Load [grams]"
colors = [:blue,:orange,:green]
lw = 2
ms = 5
msw = 0
gr()

# Beam 1, θ=0⁰
plt11 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 1, \$\\theta=0^{\\degree}\$")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label=" Experimental")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-u1_500mm[1,1], u2_500mm[1,1], -u3_500mm[1,1]], colors)
    plot!(x, σVector[1,1]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b1_th0_ref[1,:], u2_b1_th0_ref[1,:], u3_b1_th0_ref[1,:]], [u1_b1_th0_ref[2,:], u2_b1_th0_ref[2,:], u3_b1_th0_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt11)
savefig(string(absPath,"/compositeCantileverMD_b1th0.pdf"))

# Beam 1, θ=45⁰
plt12 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 1, \$\\theta=45^{\\degree}\$")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="  Experimental")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-u1_500mm[1,2], u2_500mm[1,2], -u3_500mm[1,2]], colors)
    plot!(x, σVector[1,2]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b1_th45_ref[1,:], u2_b1_th45_ref[1,:], u3_b1_th45_ref[1,:]], [u1_b1_th45_ref[2,:], u2_b1_th45_ref[2,:], u3_b1_th45_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt12)
savefig(string(absPath,"/compositeCantileverMD_b1th45.pdf"))

# Beam 2, θ=0⁰
plt21 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 2, \$\\theta=0^{\\degree}\$")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="  Experimental")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-u1_500mm[2,1], u2_500mm[2,1], -u3_500mm[2,1]], colors)
    plot!(x, σVector[2,1]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b2_th0_ref[1,:], u2_b2_th0_ref[1,:], u3_b2_th0_ref[1,:]], [u1_b2_th0_ref[2,:], u2_b2_th0_ref[2,:], u3_b2_th0_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt21)
savefig(string(absPath,"/compositeCantileverMD_b2th0.pdf"))

# Beam 2, θ=45⁰
plt22 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 2, \$\\theta=45^{\\degree}\$")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="  Experimental")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-u1_500mm[2,2], u2_500mm[2,2], -u3_500mm[2,2]], colors)
    plot!(x, σVector[2,2]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b2_th45_ref[1,:], u2_b2_th45_ref[1,:], u3_b2_th45_ref[1,:]], [u1_b2_th45_ref[2,:], u2_b2_th45_ref[2,:], u3_b2_th45_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
plot!(legend=:bottomright)
display(plt22)
savefig(string(absPath,"/compositeCantileverMD_b2th45.pdf"))

# Beam 3, θ=0⁰
plt31 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 3, \$\\theta=0^{\\degree}\$")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="  Experimental")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-u1_500mm[3,1], u2_500mm[3,1], -u3_500mm[3,1]], colors)
    plot!(x, σVector[3,1]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b3_th0_ref[1,:], u2_b3_th0_ref[1,:], u3_b3_th0_ref[1,:]], [u1_b3_th0_ref[2,:], u2_b3_th0_ref[2,:], u3_b3_th0_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt31)
savefig(string(absPath,"/compositeCantileverMD_b3th0.pdf"))

# Beam 3, θ=45⁰
plt32 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 3, \$\\theta=45^{\\degree}\$")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="  Experimental")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-u1_500mm[3,2], u2_500mm[3,2], -u3_500mm[3,2]], colors)
    plot!(x, σVector[3,2]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b3_th45_ref[1,:], u2_b3_th45_ref[1,:], u3_b3_th45_ref[1,:]], [u1_b3_th45_ref[2,:], u2_b3_th45_ref[2,:], u3_b3_th45_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt32)
savefig(string(absPath,"/compositeCantileverMD_b3th45.pdf"))

println("Finished compositeCantileverMDPlotGenerator.jl")