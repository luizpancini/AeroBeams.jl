using AeroBeams, Plots

pTest = rotation_between_WM([0;0;0],[0;-4.1;0])

x = -16:0.01:16

p1 = [0; 0; 0]
p2 = [[0; x[i]; 0] for i in eachindex(x)]

pH = Vector{Vector{Float64}}(undef,length(x))
ϕ = Vector{Float64}(undef,length(x))
for i in eachindex(x)
    pH[i] = rotation_between_WM(p1,p2[i])
    ϕ[i] = rotation_angle(pH[i])
end

plt1 = plot(xlabel="p2", ylabel="pH", xlims=[x[1],x[end]], xticks=x[1]:4:x[end])
plot!(x, [pH[i][2] for i in eachindex(x)], lw=2, c=:black, label=false)
display(plt1)

plt2 = plot(xlabel="p2", ylabel="ϕ [deg]", xlims=[x[1],x[end]], xticks=x[1]:4:x[end], yticks=-180:90:180)
plot!(x, ϕ*180/pi, lw=2, c=:black, label=false)
display(plt2)