using LinearInterpolations

# Define the grid points and values
x = [1.0, 2.0, 3.0]
y = [4.0, 5.0, 6.0]
z = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

# Create the interpolation object
interp_func = Interpolate((x,y), z)

# Use the interpolation function
xi, yi = 2.5, 5.5
zi = interp_func([xi, yi])
println("Interpolated value at ($xi, $yi): $zi")
