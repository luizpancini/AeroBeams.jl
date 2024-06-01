# Define the structures
struct Spring
    node::Int
end

struct Beam
    springs::Vector{Spring}
end

# Sample data: Vector of Beams with their Springs
beams = [
    Beam([Spring(1), Spring(2)]),
    Beam([Spring(3), Spring(4)]),
    Beam([Spring(5)])
]

# Extracting all nodes from all springs of all beams
all_nodes = [spring.node for beam in beams for spring in beam.springs]

# Print the result
println(all_nodes)
