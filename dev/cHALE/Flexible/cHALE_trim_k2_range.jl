using AeroBeams, JLD2, Dierckx

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 1

# Altitude
h = 20e3

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Discretization
nElemWing = 80
nElemTailBoom = 1
nElemHorzStabilizer = 2
nElemVertStabilizer = 1

# System solver
relaxFactor = 0.5
maxIter = 1000
σ0 = 1.0
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,pseudoInverseMethod=:dampedLeastSquares,displayStatus=false)

# Set bending curvature and airspeed ranges
k2Range = range(-0.015, 0.045, 5)
URange = unique(sort(vcat(20:0.5:50,42.3,42.4)))

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(k2Range),length(URange))
highestConvUindex = Array{Int64}(undef,length(k2Range))

trimAoA = fill(NaN, length(k2Range), length(URange))
trimThrust = fill(NaN, length(k2Range), length(URange))
trimδ = fill(NaN, length(k2Range), length(URange))
trimEnduranceFactor = fill(NaN, length(k2Range), length(URange))

x1_0 = Array{Vector{Float64}}(undef,length(k2Range))
x2_0 = Array{Vector{Float64}}(undef,length(k2Range))
x3_0 = Array{Vector{Float64}}(undef,length(k2Range))
x1_n = Array{Vector{Float64}}(undef, length(k2Range))
x1_e = Array{Vector{Float64}}(undef, length(k2Range))
x1_def = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
x2_def = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
x3_def = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
twist = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
cl = Array{Vector{Float64}}(undef,length(k2Range),length(URange))

# ELement ranges
elemRangeRightWing = 1 + div(nElemWing,2) : nElemWing

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for k2 = $k2, U = $U m/s")
        # Model for trim problem
        model,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
        # Set initial guess solution as previous known solution
        x0Trim = j == 1 ? zeros(0) : trimProblem[i,j-1].x
        # Create and solve trim problem
        trimProblem[i,j] = create_TrimProblem(model=model,systemSolver=NR,x0=x0Trim)
        solve!(trimProblem[i,j])
        # Skip if unconverged
        if !trimProblem[i,j].systemSolver.convergedFinalSolution
            highestConvUindex[i] = j-1
            break
        else
            highestConvUindex[i] = j
        end
        # Extract trim variables
        trimAoA[i,j] = trimProblem[i,j].aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
        trimThrust[i,j] = stabilizersAero ? trimProblem[i,j].x[end-1]*trimProblem[i,j].model.forceScaling : trimProblem[i,j].x[end]*trimProblem[i,j].model.forceScaling
        trimδ[i,j] = stabilizersAero ? trimProblem[i,j].x[end] : 0
        println("Trim AoA = $(trimAoA[i,j]*180/π), trim thrust = $(trimThrust[i,j]), trim δ = $(trimδ[i,j]*180/π)")
        # Trim cl^1.5/cd (endurance factor)
        lift = model.mass*model.atmosphere.g - trimThrust[i,j]*sin(trimAoA[i,j])
        drag = trimThrust[i,j]*cos(trimAoA[i,j])
        qS = 1/2*model.atmosphere.ρ*U^2*(32*1)
        trimEnduranceFactor[i,j] = lift^1.5/drag*sqrt(1/qS)
        # Undeformed jig-shape properties
        if j == 1
            # Undeformed nodal positions of right wing
            x1_0[i] = vcat([vcat(model.elements[e].r_n1[1],model.elements[e].r_n2[1]) for e in elemRangeRightWing]...)
            x2_0[i] = vcat([vcat(model.elements[e].r_n1[2],model.elements[e].r_n2[2]) for e in elemRangeRightWing]...)
            x3_0[i] = vcat([vcat(model.elements[e].r_n1[3],model.elements[e].r_n2[3]) for e in elemRangeRightWing]...)
            # Nodal arclength positions
            x1_n[i] = vcat([vcat(model.elements[e].x1_n1,model.elements[e].x1_n2) for e in elemRangeRightWing]...)
            # Elemental arclength positions
            x1_e[i] = vcat([model.elements[e].x1 for e in elemRangeRightWing]...)
        end
        # Displacements over span of right wing
        u1_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].u_n1[1],trimProblem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in elemRangeRightWing]...)
        u2_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].u_n1[2],trimProblem[i,j].nodalStatesOverσ[end][e].u_n2[2]) for e in elemRangeRightWing]...)
        u3_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].u_n1[3],trimProblem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in elemRangeRightWing]...)
        u1_of_x1 .-= u1_of_x1[1]
        u2_of_x1 .-= u2_of_x1[1]
        u3_of_x1 .-= u3_of_x1[1]
        # Deformed nodal positions of right wing
        x1_def[i,j] = x1_0[i] .+ u1_of_x1
        x2_def[i,j] = x2_0[i] .+ u2_of_x1
        x3_def[i,j] = x3_0[i] .+ u3_of_x1
        # Angle of twist
        p1_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].p_n1[1],trimProblem[i,j].nodalStatesOverσ[end][e].p_n2[1]) for e in elemRangeRightWing]...)
        p2_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].p_n1[2],trimProblem[i,j].nodalStatesOverσ[end][e].p_n2[2]) for e in elemRangeRightWing]...)
        p3_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].p_n1[3],trimProblem[i,j].nodalStatesOverσ[end][e].p_n2[3]) for e in elemRangeRightWing]...)
        twist[i,j] = [asind((first(rotation_tensor_WM([p1_of_x1[k],p2_of_x1[k],p3_of_x1[k]]))*AeroBeams.a2)[3]) for k in eachindex(p1_of_x1)]
        twist[i,j] .-= twist[i,j][1] # discount root angle (rigid-body rotation)
        # Lift distribution over the wing
        α = vcat([trimProblem[i,j].aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in elemRangeRightWing]...)
        cn = vcat([trimProblem[i,j].aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in elemRangeRightWing]...)
        ct = vcat([trimProblem[i,j].aeroVariablesOverσ[end][e].aeroCoefficients.ct for e in elemRangeRightWing]...)
        cl[i,j] = @. cn*cos(α) + ct*sin(α)
    end
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALE_trim_k2_range"
absPath = string(pwd(),relPath)
mkpath(absPath)
relPathData = "/dev/cHALE/Flexible/outputs/data/cHALE_trim_k2_range/"
absPathData = string(pwd(),relPathData)
mkpath(absPathData)

# Save solution for k2 = 0.045 around the discontinuity
if λ == 1 && 0.045 in k2Range
    ik20045 = findfirst(x->x≈0.045,k2Range)
    jU424 = findfirst(x->x≈42.4,URange)
    trimProblem_k2_0045_U424 = trimProblem[ik20045,jU424]
    @save absPathData*string("cHALE_trim_lambda",λ,"_k20045_U424.jld2") trimProblem_k2_0045_U424
end

# Plot configurations
colors = palette([:royalblue, :blueviolet, :deeppink, :darkorange, :gold])
colorsGrad = cgrad(colors, 1e2, categorical=false)
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
msw = 0
mshape = [:circle, :star, :utriangle, :pentagon]
labels = ["\$k_2 = $(k2) \$" for k2 in k2Range]
L = 16
gr()

# Plot du3/dx1
for (i,k2) in enumerate(k2Range)
    plt_du3_x1 = plot(xlabel="Normalized arclength", ylabel="\$\\partial u_3 / \\partial x_1\$", title=labels[i], xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
    for (j,U) in enumerate(URange)
        x = x1_n[i][1:2:end]
        y = x3_def[i,j][1:2:end]
        spl = Spline1D(x, y; k=3)
        dy  = derivative(spl, x)
        plot!(x/L, dy, lw=lw, lz=U, c=colorsGrad, label=false, colorbar_title="Airspeed [m/s]")
    end
    display(plt_du3_x1)
    savefig(string(absPath,"/cHALE_trim_k2_range_du3dx1_lambda",λ,"_k2_",k2,".pdf"))
end

# c_l and lift distributions over airspeed for selected k2
for (i,k2) in enumerate(k2Range)
    plt_cl_dist = plot(xlabel="Normalized arclength", ylabel="\$c_l\$", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
    for (j,U) in enumerate(URange)
        plot!(x1_e[i]/L, cl[i,j], lz=U, c=colorsGrad, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
    end
    display(plt_cl_dist)
    savefig(string(absPath,"/cHALE_trim_k2_range_cl_dist_lambda",λ,"_k2_",k2,".pdf"))

    plt_L_dist = plot(xlabel="Normalized arclength", ylabel="Lift/span [N/m]", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
    for (j,U) in enumerate(URange)
        plot!(x1_e[i]/L, 1/2*trimProblem[i,j].model.atmosphere.ρ*U^2*trimProblem[i,j].model.beams[2].aeroSurface.c*cl[i,j], lz=U, c=colorsGrad, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
    end
    display(plt_L_dist)
    savefig(string(absPath,"/cHALE_trim_k2_range_L_dist_lambda",λ,"_k2_",k2,".pdf"))
end

# Normalized deformed span in x1-x3 plane at airspeed extrema
plt_disp = plot(xlabel="Normalized spanwise direction", ylabel="Normalized vertical direction", xlims=[0,1], ylims=[0,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=(0.4,0.9))
plot!([NaN], [NaN], ls=:solid, c=:black, lw=lw, label="\$U_{\\infty} = $(URange[1])\\;\\mathrm{m/s}\$")
plot!([NaN], [NaN], ls=:dash, c=:black, lw=lw, label="\$U_{\\infty} = $(URange[end])\\;\\mathrm{m/s}\$")
for (i,k2) in enumerate(k2Range)
    plot!(x1_def[i,1]/L, x3_def[i,1]/L, ls=:solid, c=colors[i], lw=lw, label=labels[i])
    plot!(x1_def[i,end]/L, x3_def[i,end]/L, ls=:dash, c=colors[i], lw=lw, label=false)
end
display(plt_disp)
savefig(string(absPath,"/cHALE_trim_k2_range_disp_lambda",λ,".pdf"))

# Normalized deformed span in x1-x2 plane at airspeed extrema
plt_dispIP = plot(xlabel="Normalized spanwise direction", ylabel="Normalized horizontal direction", xlims=[0,1], ylims=[-1,0], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:best)
plot!([NaN], [NaN], ls=:solid, c=:black, lw=lw, label="\$U_{\\infty} = $(URange[1])\\;\\mathrm{m/s}\$")
plot!([NaN], [NaN], ls=:dash, c=:black, lw=lw, label="\$U_{\\infty} = $(URange[end])\\;\\mathrm{m/s}\$")
for (i,k2) in enumerate(k2Range)
    plot!(x1_def[i,1]/L, x2_def[i,1]/L, ls=:solid, c=colors[i], lw=lw, label=labels[i])
    plot!(x1_def[i,end]/L, x2_def[i,end]/L, ls=:dash, c=colors[i], lw=lw, label=false)
end
display(plt_dispIP)
savefig(string(absPath,"/cHALE_trim_k2_range_dispIP_lambda",λ,".pdf"))

# Angle of twist at airspeed extrema
plt_twist = plot(xlabel="Normalized arclength", ylabel="Angle of twist [deg]", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(x1_n[i]/L, twist[i,1], c=colors[i], lw=lw, ls=:solid, label=false)
    plot!(x1_n[i]/L, twist[i,end], c=colors[i], lw=lw, ls=:dash, label=false)
end
display(plt_twist)
savefig(string(absPath,string("/cHALE_trim_k2_range_twist_lambda",λ,".pdf")))

# Trim root angle of attack
plt_trimAoA = plot(xlabel="Airspeed [m/s]", ylabel="Trim root angle of attack [deg]", xlims=extrema(URange), ylims=[-5,20], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topright)
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimAoA[i,:]*180/π, c=colors[i], lw=lw, label=false)
end
display(plt_trimAoA)
savefig(string(absPath,"/cHALE_trim_k2_range_trimAoA_lambda",λ,".pdf"))

# Trim thrust
plt_trimT = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=extrema(URange), ylims=[0,100], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimThrust[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_trimT)
savefig(string(absPath,"/cHALE_trim_k2_range_trimThrust_lambda",λ,".pdf"))

# Trim elevator deflection
plt_trimDelta = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=extrema(URange), ylims=[-50,10], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimδ[i,:]*180/π, c=colors[i], lw=lw, label=false)
end
display(plt_trimDelta)
savefig(string(absPath,"/cHALE_trim_k2_range_trimDelta_lambda",λ,".pdf"))

# Endurance factor (cl^1.5/cd)
plt_enduranceFactor = plot(xlabel="Airspeed [m/s]", ylabel="Trim \$c_L^{1.5}/c_D\$", xlims=extrema(URange), ylims=[0,40], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimEnduranceFactor[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_enduranceFactor)
savefig(string(absPath,"/cHALE_trim_k2_range_enduranceFactor_lambda",λ,".pdf"))

println("Finished cHALE_trim_k2_range.jl")