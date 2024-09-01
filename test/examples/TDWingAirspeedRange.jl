using AeroBeams, LinearAlgebra, DelimitedFiles

# Discretization
nElem = 20

# Loop configurations
θRange = π/180*[1.0; 2.2]
URange = LinRange(0,40,41)
freqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
tip_u3 = Array{Float64}(undef,length(θRange),length(URange))
tip_twist = Array{Float64}(undef,length(θRange),length(URange))
for (i,θ) in enumerate(θRange)
    # Model
    TDWingAirspeedRange = create_TDWing(θ=θ,nElem=nElem)
    # Loop airspeeds
    for (j,U) in enumerate(URange)
        display("Solving for root angle = $(θ*180/π) deg, U = $(round(Int,U)) m/s")
        # Update velocity of basis A 
        set_motion_basis_A!(model=TDWingAirspeedRange,v_A=[0;U;0])
        # Create and solve problem
        global problem = create_EigenProblem(model=TDWingAirspeedRange,nModes=4,frequencyFilterLimits=[0.1,Inf64],normalizeModeShapes=true)
        solve!(problem)
        # Get outputs
        freqs[i,j] = problem.frequenciesOscillatory/(2π)
        tip_u3[i,j] = problem.nodalStatesOverσ[end][nElem].u_n2[3]
        tip_p = problem.nodalStatesOverσ[end][nElem].p_n2
        R,_ = rotation_tensor_WM(tip_p)
        Δ = R*[0; 1; 0]
        tip_twist[i,j] = asind(Δ[3])
    end
end

# Load reference solutions
u3_1deg_exp = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "u3_1deg_exp.txt"))
u3_1deg_num = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "u3_1deg_num.txt"))
u3_2_2deg_exp = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "u3_2_2deg_exp.txt"))
u3_2_2deg_num = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "u3_2_2deg_num.txt"))
th_1deg_exp = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "th_1deg_exp.txt"))
th_1deg_num = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "th_1deg_num.txt"))
th_2_2deg_exp = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "th_2_2deg_exp.txt"))
th_2_2deg_num = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "th_2_2deg_num.txt"))
freqs_ref = readdlm(joinpath(dirname(@__DIR__), "referenceData", "TDWing", "freqs.txt"))

println("Finished TDWingAirspeedRange.jl")