using DelimitedFiles

# Nested results use one row per parameter combination, in column-major order
test_data_matrix(data::Number) = reshape([data], 1, 1)
test_data_matrix(data::AbstractArray{<:Number}) = ndims(data) <= 2 ? reshape(data, size(data, 1), size(data, 2)) : reshape(vec(data), :, 1)

_flatten_test_data(data::Number) = [Float64(data)]
_flatten_test_data(data::AbstractArray{<:Number}) = Float64.(vec(data))
function _flatten_test_data(data::AbstractArray{<:AbstractArray})
    flattened = Float64[]
    for nestedData in data
        append!(flattened, _flatten_test_data(nestedData))
    end
    return flattened
end

function test_data_matrix(data::AbstractArray{<:AbstractArray})
    flattened = [_flatten_test_data(row) for row in vec(data)]
    isempty(flattened) && return reshape(Float64[], 0, 0)
    maxLength = maximum(length.(flattened))
    matrix = fill(NaN, length(flattened), maxLength)
    for i in eachindex(flattened)
        matrix[i, 1:length(flattened[i])] .= flattened[i]
    end
    return matrix
end

# Helper function to match test data against reference data
function matches_test_data(data, example, filename; atol=SELFatol, rtol=SELFrtol)
    path = joinpath(@__DIR__, "newTestDataGenerators", example, filename)
    isfile(path) || throw(ArgumentError("Reference data file not found: $path"))
    # DelimitedFiles cannot read an empty file
    if filesize(path) == 0
        return isempty(data)
    end
    reference = readdlm(path, Float64)
    actual = test_data_matrix(data)
    return size(actual) == size(reference) && isapprox(actual, reference; atol, rtol, nans=true)
end

# Include helper for example scripts, with optional skipping of plot/visualization sections
function include_example(path::AbstractString)
    if !runExamplePlots
        sourcePath = joinpath(@__DIR__, path)
        sourceCode = read(sourcePath, String)
        markerMatch = findfirst(r"(?m)^##\s*(Plots|Plot|Visualization)\b", sourceCode)
        if markerMatch !== nothing
            sourceCode = sourceCode[1:prevind(sourceCode, first(markerMatch))]
            return Base.include_string(@__MODULE__, sourceCode, sourcePath)
        end
    end
    return include(path)
end