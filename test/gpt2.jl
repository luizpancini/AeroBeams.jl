using Plots

x1 = collect(1:0.5:3)
x2 = collect(1.5:0.5:3.5)
y1 = collect(3:-0.5:1)
y2 = collect(3:-0.5:1)

plt = plot(Shape(vcat(x1,reverse(x2)),vcat(y1,reverse(y2))), fillcolor = plot_color(:red, 0.25), lw=2, label="Flutter region")
display(plt)