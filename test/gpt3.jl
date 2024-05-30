using Plots, Statistics

gr(size=(600,600))

function plot_iso3d(xs, ys, zs; lw=3, lc=:red, title="Isometric 3D plot",label=false, camera=(45,30))
    # condition data for nearly isometric 3D plot 
    x12, y12, z12 = extrema(xs), extrema(ys), extrema(zs)
    d = maximum([diff([x12...]),diff([y12...]),diff([z12...])])[1] / 2
    xm, ym, zm = mean(x12),  mean(y12),  mean(z12) 

    # plot data
    p = plot(; xlabel="x",ylabel="y",zlabel="z", aspect_ratio=:equal, grid=:true)
    plot!(xlims=(xm-d,xm+d), ylims=(ym-d,ym+d), zlims=(zm-d,zm+d))
    plot!(;camera=camera)    #(azimuth,elevation) ???
    plot!(xs, ys, zs, title=title,lw=lw,lc=lc,label=label)
    plot!(xs, ys, zlims(p)[1] .+ 0*zs, lw=1, lc=:lightgray, label=false)
    plot!(xs, ylims(p)[2]  .+ 0*ys, zs, lw=1, lc=:lightgray, label=false)
    plot!(xlims(p)[1]  .+ 0*xs, ys, zs, lw=1, lc=:lightgray, label=false)
end

#input data
N = 100
xs = LinRange(0,100,N)
ys = LinRange(0,20,N) .+ 10*sin.(xs)
zs = LinRange(-10,50,N) .+ 15*sin.(xs)
plot_iso3d(xs, ys, zs)