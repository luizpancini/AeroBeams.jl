using AeroBeams

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
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# System solver
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)

# Set bending curvature and airspeed ranges
k2Range = range(-0.015, 0.045, 5)
URange = 20:0.5:40

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(k2Range),length(URange))
highestConvUindex = Array{Int64}(undef,length(k2Range))

trimAoA = fill(NaN, length(k2Range), length(URange))
trimThrust = fill(NaN, length(k2Range), length(URange))
trimδ = fill(NaN, length(k2Range), length(URange))
trimEnduranceFactor = fill(NaN, length(k2Range), length(URange))

x1_0 = Array{Vector{Float64}}(undef,length(k2Range))
x3_0 = Array{Vector{Float64}}(undef,length(k2Range))
x1_n = Array{Vector{Float64}}(undef, length(k2Range))
x1_def = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
x3_def = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
twist = Array{Vector{Float64}}(undef,length(k2Range),length(URange))

# ELement ranges
elemRangeRightWing = 1 + div(nElemWing,2) : nElemWing

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for k2 = $k2, U = $U m/s")
        # Model for trim problem
        model,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
        # Update model
        model.skipValidationMotionBasisA = true
        update_model!(model)
        # Set initial guess solution as previous known solution
        x0Trim = j == 1 ? zeros(0) : trimProblem[i,j-1].x
        # Create and trim problem
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
            x3_0[i] = vcat([vcat(model.elements[e].r_n1[3],model.elements[e].r_n2[3]) for e in elemRangeRightWing]...)
            # Nodal arclength positions
            x1_n[i] = vcat([vcat(model.elements[e].x1_n1,model.elements[e].x1_n2) for e in elemRangeRightWing]...)
        end
        # Displacements over span of right wing
        u1_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].u_n1[1],trimProblem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in elemRangeRightWing]...)
        u3_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].u_n1[3],trimProblem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in elemRangeRightWing]...)
        u1_of_x1 .-= u1_of_x1[1]
        u3_of_x1 .-= u3_of_x1[1]
        # Deformed nodal positions of right wing
        x1_def[i,j] = x1_0[i] .+ u1_of_x1
        x3_def[i,j] = x3_0[i] .+ u3_of_x1
        # Angle of twist
        p1_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].p_n1[1],trimProblem[i,j].nodalStatesOverσ[end][e].p_n2[1]) for e in elemRangeRightWing]...)
        p2_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].p_n1[2],trimProblem[i,j].nodalStatesOverσ[end][e].p_n2[2]) for e in elemRangeRightWing]...)
        p3_of_x1 = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].p_n1[3],trimProblem[i,j].nodalStatesOverσ[end][e].p_n2[3]) for e in elemRangeRightWing]...)
        twist[i,j] = [asind((first(rotation_tensor_WM([p1_of_x1[k],p2_of_x1[k],p3_of_x1[k]]))*AeroBeams.a2)[3]) for k in eachindex(p1_of_x1)]
        twist[i,j] .-= twist[i,j][1] # discount root angle (rigid-body rotation)
    end
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALE_trim_k2_range"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(k2Range), categorical=true)
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

# Find indices for U = 20 and U = 40
iU20 = findfirst(x->x≈20,URange)
iU40 = findfirst(x->x≈40,URange)

# Normalized deformed span at U = 20 and U = 40
plt_disp = plot(xlabel="Normalized spanwise direction", ylabel="Normalized vertical direction", xlims=[0,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topleft)
plot!([NaN], [NaN], ls=:solid, c=:black, lw=lw, label=string("\$U_{\\infty} = 20 \$ m/s"))
plot!([NaN], [NaN], ls=:dashdot, c=:black, lw=lw, label=string("\$U_{\\infty} = 40 \$ m/s"))
for (i,k2) in enumerate(k2Range)
    plot!(x1_def[i,iU20]/L, x3_def[i,iU20]/L, ls=:solid, c=colors[i], lw=lw, label=labels[i])
    plot!(x1_def[i,iU40]/L, x3_def[i,iU40]/L, ls=:dashdot, c=colors[i], lw=lw, label=false)
end
display(plt_disp)
savefig(string(absPath,"/cHALE_trim_k2_range_disp.pdf"))

# Angle of twist at U = 20 and U = 40
plt_twist = plot(xlabel="Normalized arclength", ylabel="Angle of twist [deg]", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(x1_n[i]/L, twist[i,iU20], c=colors[i], lw=lw, ls=:solid, label=false)
    plot!(x1_n[i]/L, twist[i,iU40], c=colors[i], lw=lw, ls=:dashdot, label=false)
end
display(plt_twist)
savefig(string(absPath,string("/cHALE_trim_k2_range_twist.pdf")))

# Trim root angle of attack
plt_trimAoA = plot(xlabel="Airspeed [m/s]", ylabel="Trim root angle of attack [deg]", xlims=[20,40], ylims=[-5,20], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:topright)
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimAoA[i,:]*180/π, c=colors[i], lw=lw, label=labels[i])
end
display(plt_trimAoA)
savefig(string(absPath,"/cHALE_trim_k2_range_trimAoA.pdf"))

# Trim thrust
plt_trimT = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[20,40], ylims=[0,50], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimThrust[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_trimT)
savefig(string(absPath,"/cHALE_trim_k2_range_trimThrust.pdf"))

# Trim elevator deflection
plt_trimDelta = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[20,40], ylims=[-50,10], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimδ[i,:]*180/π, c=colors[i], lw=lw, label=false)
end
display(plt_trimDelta)
savefig(string(absPath,"/cHALE_trim_k2_range_trimDelta.pdf"))

# Endurance factor (cl^1.5/cd)
plt_enduranceFactor = plot(xlabel="Airspeed [m/s]", ylabel="Trim \$c_L^{1.5}/c_D\$", xlims=[20,40], ylims=[0,40], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimEnduranceFactor[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_enduranceFactor)
savefig(string(absPath,"/cHALE_trim_k2_range_enduranceFactor.pdf"))

println("Finished cHALE_trim_k2_range.jl")