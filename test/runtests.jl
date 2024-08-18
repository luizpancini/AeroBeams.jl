using Test

@testset "Cantilever with distributed follower force static analysis" begin
    include("distributedLoadCantilever.jl")
    @test isapprox(-tip_u1[end]/L,1.4299, atol=1e-4)
    @test isapprox(tip_u3[end]/L,0.3725, atol=1e-4)
    @test isapprox(-tip_angle[end]/π,1.0321, atol=1e-4)
end

@testset "Cantilever with tip moment static analysis" begin
    include("tipMomentCantilever.jl")
    @test isapprox(tip_u1[end]/L,-1, atol=1e-3)
    @test isapprox(tip_u3[end]/L,0, atol=1e-6)
    @test isapprox(tip_angle[end]/π,2, atol=1e-4)
end

@testset "2-story frame eigen analysis" begin
    include("twoStoryFrame.jl")
    @test isapprox(freqs[[1,4]]/(2π),[11.8; 34.1], rtol=1e-2)
end

@testset "Straight rotor eigen analysis" begin
    include("unsweptTipRotor.jl")
    @test isapprox(numFreqs[end],[93.1083; 275.2196; 288.9760; 392.5491; 659.2405; 1159.1923], atol=1e-2)
end