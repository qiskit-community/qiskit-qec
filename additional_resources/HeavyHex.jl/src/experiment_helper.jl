"""
module ExperimentHelper

Contains the bit_order_functions for each job. These functions, literally, reorder the bit readout 
    based on factors in the original experiment. 
"""
module ExperimentHelper

export bit_order_function, bit_order_function_child, bit_order_function_x

function bit_order_function(r)
    return i -> bit_order_function_child(r, i)
end

function bit_order_function_child(r, i) #  test=[0, "normal bof"])
    n = 4 * (r - 1) + 2 * r
    i == n - 1 && return n - 1
    i == n - 2 && return n - r - 1
    q, _r = divrem(i, 6)
    val = (r - 1) * _r + q
    _r == 5 && return val + 1
    return val
end

function bit_order_function_x(r::Int)
    return i -> bit_order_function_x_child(r, i)
end

function bit_order_function_x_child(r::Int, i::Int) # test=[0, "normal bofX"]):
    q, _r = divrem(i, 6)
    val = r * _r + q
    _r == 5 && return val - 1
    return val
end

end # module ExperimentHelper
