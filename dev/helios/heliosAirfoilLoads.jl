using AeroBeams, ForwardDiff, DelimitedFiles

# Frames and respective data
frames = ["02_000_101", "02_010_151", "02_010_351"]
airfoil_frames, a₀_frames, a₁_frames, b_frames, k_frames, Ma_frames, U_frames = 
map(GU_frames_loader, frames) |> x -> [getindex.(x, i) for i in 1:7]

# Circulatory indicial function
circulatoryIndicialFunction = "Wagner"

# Aerodynamic solvers
aeroSolvers = [Indicial(circulatoryIndicialFunction=circulatoryIndicialFunction); BLi(circulatoryIndicialFunction=circulatoryIndicialFunction)]

# Rigid beam properties
L = 0.1
ρA,ρI = 1e-1, 1e-1
nElem = 1
∞ = 1e10

# Initialize outputs
problem = Array{DynamicProblem}(undef,length(frames),length(aeroSolvers))
τ = Array{Float64}(undef,length(frames))
t = Array{Vector{Float64}}(undef,length(frames),length(aeroSolvers))
α = Array{Vector{Float64}}(undef,length(frames),length(aeroSolvers))
cn = Array{Vector{Float64}}(undef,length(frames),length(aeroSolvers))
cm = Array{Vector{Float64}}(undef,length(frames),length(aeroSolvers))
ct = Array{Vector{Float64}}(undef,length(frames),length(aeroSolvers))

# Loop frames
for (i,frame,airfoil,a₀,a₁,b,k,Ma,U) in zip(1:length(frames),frames,airfoil_frames,a₀_frames,a₁_frames,b_frames,k_frames,Ma_frames,U_frames)
    # Create offset for middle pitch angle
    Δ = airfoil.attachedFlowParameters.α₀N
    a₀ += Δ
    # Loop aerodynamic solvers
    for (j,aeroSolver) in enumerate(aeroSolvers)
        # Display progress
        display("Solving for frame $frame, $(aeroSolver.name) solver")
        # Pitch profile
        ω = k*U/b
        τ[i] = 2π/ω
        t₀ = -τ[i]/4
        θ = t -> a₁*(1+sin(ω*(t+t₀)))*((t/τ[i])^10/(1+(t/τ[i])^10))
        p = t -> 4*tan(θ(t)/4)
        pdot = t -> ForwardDiff.derivative(p,t)
        # Update airfoil parameters
        update_Airfoil_params!(airfoil,Ma=Ma,U=U,b=b)
        # Aerodynamic surface
        surf = create_AeroSurface(solver=aeroSolver,airfoil=airfoil,c=2*b,normSparPos=0.25,updateAirfoilParameters=false)
        # Wing
        wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=ρA,ρIy=ρI,ρIz=ρI)],rotationParametrization="E321",p0=[0;0;a₀-a₁],aeroSurface=surf,pdot0_of_x1=x1->[pdot(0); 0.0; 0.0])
        # BCs
        driver = create_BC(name="driver",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,t -> p(t),0,0])
        journal = create_BC(name="journal",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])
        # Model
        model = create_Model(name="pitchingAirfoil",beams=[wing],BCs=[driver,journal],v_A=[0;U;0])
        # Time variables
        nCycles = 6
        Δt = τ[i]/1000
        tf = nCycles*τ[i]
        # Initial velocities update options
        initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2, tol=1e-12, displayProgress=true, relaxFactor=0.5, Δt=Δt/1e3)
        # Create and solve dynamic problem
        problem[i,j] = create_DynamicProblem(model=model,finalTime=tf,Δt=Δt,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
        solve!(problem[i,j])
        # Get outputs
        t[i,j] = problem[i,j].savedTimeVector
        Nt = length(t[i,j])
        α[i,j] = [problem[i,j].aeroVariablesOverTime[k][1].flowAnglesAndRates.α-Δ for k in 1:Nt]
        cn[i,j] = [problem[i,j].aeroVariablesOverTime[k][1].aeroCoefficients.cn for k in 1:Nt]
        cm[i,j] = [problem[i,j].aeroVariablesOverTime[k][1].aeroCoefficients.cm for k in 1:Nt]
        ct[i,j] = [problem[i,j].aeroVariablesOverTime[k][1].aeroCoefficients.ct for k in 1:Nt]
    end
end

# Load reference data
cnRef = Array{Matrix{Float64}}(undef,length(frames))
cmRef = Array{Matrix{Float64}}(undef,length(frames))
ctRef = Array{Matrix{Float64}}(undef,length(frames))
for (i,frame) in enumerate(frames)
    onlyNumbersFrameStr = replace(frame, "_" => "")
    cnRef[i] = readdlm(pkgdir(AeroBeams)*"/test/referenceData/GUframes/"*onlyNumbersFrameStr*"_cn.txt")
    cmRef[i] = readdlm(pkgdir(AeroBeams)*"/test/referenceData/GUframes/"*onlyNumbersFrameStr*"_cm.txt")
    ctRef[i] = readdlm(pkgdir(AeroBeams)*"/test/referenceData/GUframes/"*onlyNumbersFrameStr*"_ct.txt")
end

# Set paths
relPath = "/dev/helios/figures/heliosAirfoilLoads"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
using Plots, ColorSchemes
gr()
solverColors = cgrad(:rainbow, length(aeroSolvers), categorical=true)
frameColors = cgrad(:rainbow, length(frames), categorical=true)
lsSolvers = [:dash, :solid]
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 4
msw = 0
fps = 30
DPI = 300
labels = ["Attached flow" "Dynamic stall"]
miniLabels = ["AF" "DS"]
ΔtAnim = [2.5e-1; 1e-2; 5e-3]
YLIMS_cn = [(-0.4, 2.85), (0.4, 1.6), (-0.25, 2.1)]
YLIMS_cm = [(-0.2, 0.1), (-0.1, 0.1), (-0.35, 0.15)]

# Arrays of last cycle
tCycle = Array{Vector{Float64}}(undef,length(frames),length(aeroSolvers))
αCycle = Array{Vector{Float64}}(undef,length(frames),length(aeroSolvers))
cnCycle = Array{Vector{Float64}}(undef,length(frames),length(aeroSolvers))
cmCycle = Array{Vector{Float64}}(undef,length(frames),length(aeroSolvers))
ctCycle = Array{Vector{Float64}}(undef,length(frames),length(aeroSolvers))
for (i,frame,airfoil) in zip(eachindex(frames),frames,airfoil_frames)
    for (j,aeroSolver) in enumerate(aeroSolvers)
        range2plot = findfirst(x -> x > t[i,j][end]-τ[i],t[i,j]) : length(t[i,j])
        tCycle[i,j] = t[i,j][range2plot] .- t[i,j][end] .+ τ[i]
        αCycle[i,j] = α[i,j][range2plot] * 180/π
        cnCycle[i,j] = cn[i,j][range2plot]
        cmCycle[i,j] = cm[i,j][range2plot]
        ctCycle[i,j] = ct[i,j][range2plot]
    end
end

# cn vs alpha - separate figures
for (i,frame,airfoil) in zip(eachindex(frames),frames,airfoil_frames)
    plt_cn = plot(xlabel="Pitch angle [deg]", ylabel="\$c_n\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legendposition=:topleft)
    scatter!(cnRef[i][1,:], cnRef[i][2,:], c=:black, ms=ms, msw=msw, label="Green & Giuni (2017)")
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!(αCycle[i,j], cnCycle[i,j], c=solverColors[j], lw=lw, label=labels[j])
    end
    display(plt_cn)
    savefig(string(absPath,"/GU_",frame,"_cn.pdf"))
end

# cn vs alpha animations
for (i,frame,airfoil) in zip(eachindex(frames),frames,airfoil_frames)
    plt_cn_anim = plot(xlabel="Pitch angle [deg]", ylabel="\$c_n\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legendposition=:topleft, ylims=YLIMS_cn[i], dpi=DPI)
    scatter!(cnRef[i][1,:], cnRef[i][2,:], c=:black, ms=ms, msw=msw, label="Green & Giuni (2017)")
    for (j,aeroSolver) in enumerate(aeroSolvers)
        plot!([NaN], [NaN], c=solverColors[j], lw=lw, label=labels[j])
    end
    anim = @animate for (k,timeNow) in enumerate(tCycle[i,1])
        plotFrequency = ceil(Int,ΔtAnim[i]/problem[i,1].Δt)
        if k > 1 && rem(k,plotFrequency) > 0
            continue
        end
        plot(plt_cn_anim)
        title!("\$t = $(round(timeNow,digits=2))\$ s")
        for (j,aeroSolver) in enumerate(aeroSolvers)
            plot!(αCycle[i,j][1:k], cnCycle[i,j][1:k], c=solverColors[j], lw=lw, label=false)
        end
    end
    gif_handle = gif(anim, string(absPath,"/GU_",frame,"_cn.gif"), fps=fps)
    display(gif_handle)
end

# cm vs alpha - separate figures
for (i,frame,airfoil) in zip(eachindex(frames),frames,airfoil_frames)
    plt_cm = plot(xlabel="Pitch angle [deg]", ylabel="\$c_m\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
    scatter!(cmRef[i][1,:], cmRef[i][2,:], c=:black, ms=ms, msw=msw, label=false)
    for (j,aeroSolver) in enumerate(aeroSolvers)
        range2plot = findfirst(x -> x > t[i,j][end]-τ[i],t[i,j]) : length(t[i,j])
        plot!(α[i,j][range2plot] * 180/π, cm[i,j][range2plot], c=solverColors[j], lw=lw, label=false)
    end
    display(plt_cm)
    savefig(string(absPath,"/GU_",frame,"_cm.pdf"))
end

# cm vs alpha animations
for (i,frame,airfoil) in zip(eachindex(frames),frames,airfoil_frames)
    plt_cm_anim = plot(xlabel="Pitch angle [deg]", ylabel="\$c_m\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legendposition=:topleft, ylims=YLIMS_cm[i], dpi=DPI)
    scatter!(cmRef[i][1,:], cmRef[i][2,:], c=:black, ms=ms, msw=msw, label=false)
    anim = @animate for (k,timeNow) in enumerate(tCycle[i,1])
        plotFrequency = ceil(Int,ΔtAnim[i]/problem[i,1].Δt)
        if k > 1 && rem(k,plotFrequency) > 0
            continue
        end
        plot(plt_cm_anim)
        title!("\$t = $(round(timeNow,digits=2))\$ s")
        for (j,aeroSolver) in enumerate(aeroSolvers)
            plot!(αCycle[i,j][1:k], cmCycle[i,j][1:k], c=solverColors[j], lw=lw, label=false)
        end
    end
    gif_handle = gif(anim, string(absPath,"/GU_",frame,"_cm.gif"), fps=fps)
    display(gif_handle)
end


# ct vs alpha - separate figures
for (i,frame,airfoil) in zip(eachindex(frames),frames,airfoil_frames)
    plt_ct = plot(xlabel="Pitch angle [deg]", ylabel="\$c_t\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
    scatter!(ctRef[i][1,:], ctRef[i][2,:], c=:black, ms=ms, msw=msw, label=false)
    for (j,aeroSolver) in enumerate(aeroSolvers)
        range2plot = findfirst(x -> x > t[i,j][end]-τ[i],t[i,j]) : length(t[i,j])
        plot!(α[i,j][range2plot] * 180/π, ct[i,j][range2plot], c=solverColors[j], lw=lw, label=false)
    end
    display(plt_ct)
    savefig(string(absPath,"/GU_",frame,"_ct.pdf"))
end

# cn vs alpha - all-in-one figure
plt_cn_all = plot(xlabel="Pitch angle [deg]", ylabel="\$c_n\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, ylims=[-0.25,2.5], legend=:topleft)
scatter!([NaN],[NaN], c=:black, ms=ms, msw=msw, label="Green & Giuni (2017)")
plot!([NaN],[NaN], c=:black, lw=lw, ls=lsSolvers[1], label="Attached flow")
plot!([NaN],[NaN], c=:black, lw=lw, ls=lsSolvers[2], label="Dynamic stall")
for (i,frame,airfoil) in zip(eachindex(frames),frames,airfoil_frames)
    scatter!(cnRef[i][1,:], cnRef[i][2,:], c=frameColors[i], ms=ms, msw=msw, label=false)
    for (j,aeroSolver) in enumerate(aeroSolvers)
        range2plot = findfirst(x -> x > t[i,j][end]-τ[i],t[i,j]) : length(t[i,j])
        plot!(α[i,j][range2plot] * 180/π, cn[i,j][range2plot], c=frameColors[i], lw=lw, ls=lsSolvers[j], label=false)
    end
end
display(plt_cn_all)
savefig(string(absPath,"/cn.pdf"))

# cm vs alpha - all-in-one figure
plt_cm_all = plot(xlabel="Pitch angle [deg]", ylabel="\$c_m\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, ylims=[-0.4,0.1], legend=:bottomleft)
scatter!([NaN],[NaN], c=:black, ms=ms, msw=msw, label=false)
plot!([NaN],[NaN], c=:black, lw=lw, ls=lsSolvers[1], label=false)
plot!([NaN],[NaN], c=:black, lw=lw, ls=lsSolvers[2], label=false)
for (i,frame,airfoil) in zip(eachindex(frames),frames,airfoil_frames)
    scatter!(cmRef[i][1,:], cmRef[i][2,:], c=frameColors[i], ms=ms, msw=msw, label=false)
    for (j,aeroSolver) in enumerate(aeroSolvers)
        range2plot = findfirst(x -> x > t[i,j][end]-τ[i],t[i,j]) : length(t[i,j])
        plot!(α[i,j][range2plot] * 180/π, cm[i,j][range2plot], c=frameColors[i], lw=lw, ls=lsSolvers[j], label=false)
    end
end
display(plt_cm_all)
savefig(string(absPath,"/cm.pdf"))

# ct vs alpha - all-in-one figure
plt_ct_all = plot(xlabel="Pitch angle [deg]", ylabel="\$c_t\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, ylims=[-0.1,0.6], legend=:topleft)
scatter!([NaN],[NaN], c=:black, ms=ms, msw=msw, label=false)
plot!([NaN],[NaN], c=:black, lw=lw, ls=lsSolvers[1], label=false)
plot!([NaN],[NaN], c=:black, lw=lw, ls=lsSolvers[2], label=false)
for (i,frame,airfoil) in zip(eachindex(frames),frames,airfoil_frames)
    scatter!(ctRef[i][1,:], ctRef[i][2,:], c=frameColors[i], ms=ms, msw=msw, label=false)
    for (j,aeroSolver) in enumerate(aeroSolvers)
        range2plot = findfirst(x -> x > t[i,j][end]-τ[i],t[i,j]) : length(t[i,j])
        plot!(α[i,j][range2plot] * 180/π, ct[i,j][range2plot], c=frameColors[i], lw=lw, ls=lsSolvers[j], label=false)
    end
end
display(plt_ct_all)
savefig(string(absPath,"/ct.pdf"))

# Animations
for (i,frame,airfoil) in zip(eachindex(frames),frames,airfoil_frames)
    for (j,aeroSolver) in enumerate(aeroSolvers)
        anim = plot_dynamic_deformation(problem[i,j],backendSymbol=:plotlyjs,view=(88,0),plotLimits=([-0.4,0.4],[-0.4,0.4],[-0.4,0.4]),plotFrequency=ceil(Int,ΔtAnim[i]/problem[i,j].Δt),fps=fps,loadsSizeScaler=10,plotDistLoads=true,plotBCs=false,plotAxes=false,legendPos=:topright,DPI=DPI,save=true,savePath=string(relPath,"/heliosAirfoilLoads_frame",frame,"_",miniLabels[j],".gif"),displayProgress=true)
        display(anim)
    end
end

println("Finished heliosAirfoilLoads.jl")