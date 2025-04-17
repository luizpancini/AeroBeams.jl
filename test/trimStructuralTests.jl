# Trim structural problems

@testset "Trim analysis (reaction loads check) of a beam loaded at the middle" begin
    include("examples/freeBeamTrim.jl")
    # Analytical comparison
    @test trimForce ≈ trimForceAnalytical rtol=1e-4
    # Self-comparison
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "freeBeamTrim", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "freeBeamTrim", "M2.txt"))
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Trim analysis (reaction loads check) of a simply-supported beam loaded at the middle" begin
    include("examples/midLoadedBeamTrim.jl")
    # Self-comparison
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "midLoadedBeamTrim", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "midLoadedBeamTrim", "M2.txt"))
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Trim analysis (reaction loads check) of a right-angled frame" begin
    include("examples/rightAngledFrameTrim.jl")
    # Analytical comparison
    @test balanceVerticalForce ≈ balanceVerticalForceAnalytical rtol=1e-4
    @test balanceHorizontalForce ≈ balanceHorizontalForceAnalytical rtol=1e-4
    # Self-comparison
    balanceHorizontalForce_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rightAngledFrameTrim", "balanceHorizontalForce.txt"))[1]
    balanceVerticalForce_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rightAngledFrameTrim", "balanceVerticalForce.txt"))[1]
    @test balanceHorizontalForce ≈ balanceHorizontalForce_ atol=SELFatol
    @test balanceVerticalForce ≈ balanceVerticalForce_ atol=SELFatol
end

@testset "Trim analysis of a cantilever with tip force" begin
    include("examples/tipLoadedCantileverTrim.jl")
    # Analytical comparison
    @test u3[end] ≈ u3_analytical(L) rtol=1e-3
    @test F3[1] ≈ F3_analytical(0) rtol=1e-3
    @test M2[1] ≈ M2_analytical(0) rtol=1e-3
    # Self-comparison
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipLoadedCantileverTrim", "u3.txt"))
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipLoadedCantileverTrim", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipLoadedCantileverTrim", "M2.txt"))
    @test u3 ≈ u3_ atol=SELFatol
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end
GC.gc()
sleep(1)