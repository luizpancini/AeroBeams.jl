using AeroBeams

# Stiffness factor for the aircraft's wing
λ = 1e0

# Altitude
h = 20e3

# Discretization
nElem = 20

# Landing gear height
LGheight = 1

# Set bending curvature range
k2Range = range(-0.015,0.045,5)

# System solver
maxIter = 100
σ0 = 1.0
σstep = 0.5
NR = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,displayStatus=false)

# Initialize outputs
x1_0 = Array{Vector{Float64}}(undef,length(k2Range))
x3_0 = Array{Vector{Float64}}(undef,length(k2Range))
u1_of_x1 = Array{Vector{Float64}}(undef,length(k2Range))
u3_of_x1 = Array{Vector{Float64}}(undef,length(k2Range))
u1_of_x1 = Array{Vector{Float64}}(undef,length(k2Range))
u3_of_x1 = Array{Vector{Float64}}(undef,length(k2Range))
x1_def = Array{Vector{Float64}}(undef,length(k2Range))
x3_def = Array{Vector{Float64}}(undef,length(k2Range))

problem = Array{SteadyProblem}(undef,length(k2Range))

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Tip roller
    ρ = 1/k2
    θ = 16/ρ
    z = abs(k2) > 0 ? ρ*(1-cos(θ)) : 0
    dummyBeam = create_Beam(length=1,nElements=nElem,C=[AeroBeams.I6])
    roller = create_BC(name="roller",beam=dummyBeam,node=nElem+1,types=["u3A"],values=[z-LGheight])
    # Model
    global wingModel,L = create_SMW(nElem=nElem,stiffnessFactor=λ,∞=1e12,altitude=h,k2=k2,additionalBCs=[roller])
    # Create and solve steady problem
    problem[i] = create_SteadyProblem(model=wingModel,systemSolver=NR)
    solve!(problem[i])
    # Undeformed nodal positions of right wing
    x1_0[i] = vcat([vcat(wingModel.elements[e].r_n1[1],wingModel.elements[e].r_n2[1]) for e in 1:nElem]...)
    x3_0[i] = vcat([vcat(wingModel.elements[e].r_n1[3],wingModel.elements[e].r_n2[3]) for e in 1:nElem]...)
    # Displacements over span
    u1_of_x1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[1],problem[i].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
    u3_of_x1[i] = vcat([vcat(problem[i].nodalStatesOverσ[end][e].u_n1[3],problem[i].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
    u1_of_x1[i] .-= u1_of_x1[i][1]
    u3_of_x1[i] .-= u3_of_x1[i][1]
    # Deformed nodal positions
    x1_def[i] = x1_0[i] .+ u1_of_x1[i]
    x3_def[i] = x3_0[i] .+ u3_of_x1[i]
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/outputs/figures/cHALEwing_on_ground"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:rainbow], range(0, 1, length(k2Range)))
lw = 2
labels = ["\$k_2 = $(k2) \$" for k2 in k2Range]
gr()

# Normalized deformed span
plt1 = plot(xlabel="Normalized spanwise direction", ylabel="OOP direction [m]", xlims=[0,1], ylims=[-LGheight,Inf64])
for (i,k2) in enumerate(k2Range)
    plot!(x1_def[i]/L, x3_def[i], c=colors[i], lw=lw, label=labels[i])
end
display(plt1)
savefig(string(absPath,"/cHALEwing_on_ground.pdf"))

println("Finished cHALEwing_on_ground.jl")