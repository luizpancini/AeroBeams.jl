using LinearInterpolations

# Sample data: time and corresponding quantity
time = [0.0, 1.0, 2.0, 3.0, 4.0]
quantity = [0.0, 2.0, 1.0, 3.0, 4.0]

# Create the interpolation object
itp = Interpolate(time, quantity)

# Example usage
t_query = 2.5
v = t -> itp(t)
q = v(t_query)
println("The interpolated quantity at time $t_query is $q")
