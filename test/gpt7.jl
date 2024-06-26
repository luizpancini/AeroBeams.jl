# Step 1: Initialize the outer vector with a specified size
n = 5  # specify the size
outer_vector = Vector{Vector{Int}}(undef, n)

# # Step 2: Assign empty vectors to each element
# for i in 1:n
#     outer_vector[i] = Vector{Int}()
# end

# Now you can use push!() to add elements to the inner vectors
push!(outer_vector[1], 10)
push!(outer_vector[2], 20)
push!(outer_vector[2], 30)

# Display the result
println(outer_vector)
