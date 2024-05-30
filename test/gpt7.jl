# Define the parameters
trimδ = 1
Δδ = 1
t1 = 1
t2 = t1 + 1
t3 = t2 + 1

# Define the triangle pulse function
δ = t -> ifelse(
    t <= t1, 
    trimδ,  # Linearly increasing from trimδ to trimδ + Δδ
    ifelse(
        t <= t2, 
        trimδ + Δδ * ((t-t1) / (t2-t1)),  # Linearly decreasing from trimδ + Δδ to trimδ
        ifelse(
            t <= t3, 
            trimδ + Δδ - Δδ * ((t-t2) / (t3-t2)),  # Linearly decreasing from trimδ to 0
            trimδ  # Zero after t3
        )
    )
)

# Example usage
times = 0:0.1:4  # Time range for testing
values = [δ(t) for t in times]

# Print the results
plot(times,δ.(times))
