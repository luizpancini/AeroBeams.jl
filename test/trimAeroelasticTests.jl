# Aeroelastic trim problems

@testset "Trim analysis of the Blended-Wing-Body flying wing in free flight" begin
    include("examples/BWBtrim.jl")
    # Self-comparison
    trimAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBtrim", "trimAoA.txt"))
    trimThrust_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBtrim", "trimThrust.txt"))
    trimδ_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBtrim", "trimDelta.txt"))
    @test trimAoA ≈ trimAoA_ atol=SELFatol rtol=SELFrtol
    @test trimThrust ≈ trimThrust_ atol=SELFatol rtol=SELFrtol
    @test trimδ ≈ trimδ_ atol=SELFatol rtol=SELFrtol
end

@testset "Trim analysis the conventional HALE aircraft in free flight (considering aerodynamics from stabilizers and thrust)" begin
    include("examples/conventionalHALEfullTrim.jl")
    # Self-comparison
    trimAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEfullTrim", "trimAoA.txt"))
    trimThrust_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEfullTrim", "trimThrust.txt"))
    trimδ_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEfullTrim", "trimDelta.txt"))
    @test trimAoA ≈ trimAoA_ atol=SELFatol rtol=SELFrtol
    @test trimThrust ≈ trimThrust_ atol=SELFatol rtol=SELFrtol
    @test trimδ ≈ trimδ_ atol=SELFatol rtol=SELFrtol
end

@testset "Trim analysis the conventional HALE aircraft in free flight at rigid and flexible configurations (neglecting aerodynamics from stabilizers and thrust)" begin
    include("examples/conventionalHALEtrim.jl")
    # Reference comparison
    @test trimAoA[1,1] ≈ trimAoAERef[2,1] rtol=2e-2
    @test trimAoA[2,1] ≈ trimAoARRef[2,1] rtol=2e-2
    # Self-comparison
    trimAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEtrim", "trimAoA.txt"))
    @test trimAoA ≈ trimAoA_ atol=SELFatol rtol=SELFrtol
end

@testset "Trim analysis of the Helios flying-wing" begin
    include("examples/heliosTrim.jl")
    # Self-comparison
    trimAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosTrim", "trimAoA.txt"))
    trimThrust_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosTrim", "trimThrust.txt"))
    trimδ_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosTrim", "trimDelta.txt"))
    @test trimAoA ≈ trimAoA_ atol=SELFatol rtol=SELFrtol
    @test trimThrust ≈ trimThrust_ atol=SELFatol rtol=SELFrtol
    @test trimδ ≈ trimδ_ atol=SELFatol rtol=SELFrtol
end
GC.gc()
sleep(1)