using Plots

function make_batman()
    p = [(0, 0), (0.5, 0.2), (1, 0), (1, 2),  (0.3, 1.2), (0.2, 2), (0, 1.7)]
    s = [(0.2, 1), (0.4, 1), (2, 0), (0.5, -0.6), (0, 0), (0, -0.15)]
    m = [(p[i] .+ p[i + 1]) ./ 2 .+ s[i] for i in 1:length(p) - 1]

    pts = similar(m, 0)
    for (i, mi) in enumerate(m)
        append!(
            pts,
            map(BezierCurve([p[i], m[i], p[i + 1]]), range(0, 1, length = 30))
        )
    end
    x, y = Plots.unzip(Tuple.(pts))
    Shape(vcat(x, -reverse(x)), vcat(y, reverse(y)))
end

# background and limits
plt = plot(
    bg = :black,
    xlim = (0.1, 0.9),
    ylim = (0.2, 1.5),
    framestyle = :none,
    size = (400, 400),
    legend = false,
)

# create an ellipse in the sky
pts = Plots.partialcircle(0, 2π, 100, 0.1)
x, y = Plots.unzip(pts)
x = 1.5x .+ 0.7
y .+= 1.3
pts = collect(zip(x, y))

# beam
beam = Shape([(0.3, 0.0), pts[95], pts[50], (0.3, 0.0)])
plot!(beam, fillcolor = plot_color(:yellow, 0.3))