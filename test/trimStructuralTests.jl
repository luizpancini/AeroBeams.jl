@testset "Trim analysis (reaction loads check) of a beam loaded at the middle" begin
    include_example("examples/trimStructural/freeBeamTrim.jl")
    @test matches_test_data(F3, "freeBeamTrim", "F3.txt")
    @test matches_test_data(M2, "freeBeamTrim", "M2.txt")
end

@testset "Trim analysis (reaction loads check) of a simply-supported beam loaded at the middle" begin
    include_example("examples/trimStructural/midDispBeamTrim.jl")
    @test matches_test_data(F3, "midDispBeamTrim", "F3.txt")
    @test matches_test_data(M2, "midDispBeamTrim", "M2.txt")
end

@testset "Trim analysis (reaction loads check) of a right-angled frame" begin
    include_example("examples/trimStructural/rightAngledFrameTrim.jl")
    @test matches_test_data(balanceHorizontalForce, "rightAngledFrameTrim", "balanceHorizontalForce.txt")
    @test matches_test_data(balanceVerticalForce, "rightAngledFrameTrim", "balanceVerticalForce.txt")
end

@testset "Trim analysis of a cantilever with tip force" begin
    include_example("examples/trimStructural/tipLoadedCantileverTrim.jl")
    @test matches_test_data(u3, "tipLoadedCantileverTrim", "u3.txt")
    @test matches_test_data(F3, "tipLoadedCantileverTrim", "F3.txt")
    @test matches_test_data(M2, "tipLoadedCantileverTrim", "M2.txt")
end