function biased_vector(a, b, n, c)
    step = (b - a) / (1 + (c - 1) * (n - 1))

    vector = [a]
    current_value = a
    for i in 2:n
        current_value += step * c^(1 - i)
        push!(vector, current_value)
    end
    return vector
end

# Example usage
a = 0.
b = 10.
n = 5
c = 2.

biased_vec = biased_vector(a, b, n, c)
println(biased_vec)
