using ForwardDiff

# Function to convert arrays from Dual numbers to their values
function convert_to_values(arr)
    if isa(arr, AbstractArray)
        return map(convert_to_values, arr)
    elseif isa(arr, ForwardDiff.Dual)
        return ForwardDiff.value(arr)
    else
        return arr
    end
end

# Define the function to be differentiated
function my_function(x)
    return [x[1]^2 + x[2]^2, x[1] * x[2]]
end

# Original variables
x = [1.0, 2.0]

# Compute the Jacobian in a local scope
jacobian_result = ForwardDiff.jacobian(my_function, x) 

# Example usage:

# Define an array with Dual numbers
dual_array = [ForwardDiff.Dual(1.0, 1.0) ForwardDiff.Dual(2.0, 1.0); ForwardDiff.Dual(3.0, 1.0) ForwardDiff.Dual(4.0, 1.0)]

# Convert the array from Dual numbers to their values
value_array = convert_to_values(dual_array)

# Print the original Dual array and the converted value array
println("Original Dual array:")
println(dual_array)

println("Converted value array:")
println(value_array)
