@testset "Reference data layouts" begin
    matrix = [1.0 3.0 5.0; 2.0 4.0 6.0]
    nested = [matrix, matrix .+ 10]
    @test test_data_matrix(nested) == [1 2 3 4 5 6; 11 12 13 14 15 16]
    @test test_data_matrix(reshape(vcat(nested, nested), 2, 2)) ==
        repeat(test_data_matrix(nested), 2, 1)
    ragged = [[1.0 2.0; 3.0 4.0], [5.0 6.0]]
    @test isapprox(test_data_matrix(ragged), [1 3 2 4; 5 6 NaN NaN]; atol=0, rtol=0, nans=true)
    deeplyNested = [[[1.0, 2.0], [3.0, 4.0]], [[5.0], [6.0, 7.0, 8.0]]]
    @test test_data_matrix(deeplyNested) == [1 2 3 4; 5 6 7 8]
    @test test_data_matrix([1.0, 2.0]) == reshape([1.0, 2.0], 2, 1)
    @test test_data_matrix(3.0) == fill(3.0, 1, 1)

    mktempdir() do dir
        for data in (3.0, [1.0, 2.0], matrix, nested,
                deeplyNested,
                reshape(vcat(nested, nested), 2, 2),
                reshape(collect(1.0:8.0), 2, 2, 2), [1.0, NaN, Inf])
            writedlm(joinpath(dir, "data.txt"), test_data_matrix(data))
            @test matches_test_data(data, dir, "data.txt"; atol=0, rtol=0)
        end

        writedlm(joinpath(dir, "data.txt"), matrix)
        @test !matches_test_data(permutedims(matrix), dir, "data.txt")
        @test !matches_test_data(matrix .+ 1, dir, "data.txt")
        writedlm(joinpath(dir, "data.txt"), Float64[])
        @test matches_test_data(Float64[], dir, "data.txt")
        @test !matches_test_data([1.0], dir, "data.txt")
        @test_throws ArgumentError matches_test_data(matrix, dir, "missing.txt")
        writedlm(joinpath(dir, "data.txt"), ["invalid"])
        @test_throws ErrorException matches_test_data(matrix, dir, "data.txt")
        writedlm(joinpath(dir, "data.txt"), test_data_matrix(ragged))
        @test matches_test_data(ragged, dir, "data.txt"; atol=0, rtol=0)
    end
end

@testset "Optional example plots" begin
    savedRunExamplePlots = runExamplePlots
    try
        mktempdir() do dir
            path = joinpath(dir, "example.jl")
            for marker in ("## Plots", "## Plot", "## Visualization")
                write(path, "1\n$marker\n2\n")
                global runExamplePlots = false
                @test include_example(path) == 1
                global runExamplePlots = true
                @test include_example(path) == 2
            end
            write(path, "23\n")
            @test include_example(path) == 23
        end
    finally
        global runExamplePlots = savedRunExamplePlots
    end
end
