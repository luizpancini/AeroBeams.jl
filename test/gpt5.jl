using Plots
n = 10
x = collect(0:10:100)
y = rand(length(x), n)
plt = plot()
for i in 1:n
    plot!(x, y[:,i], lz=i/n, c=:rainbow, label = false, lw = 3, colorbar_title="AIRSPEED")
end
display(plt)