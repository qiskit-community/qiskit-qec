module AnalyzerHelper

# Using eval adds some performance; the value of r is a compiled constant.
# In practice, bit_order_function is a small part of the performance budget.
# bit_order_function(r) = @eval i -> bit_order_function($r, i)

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

end #module AnalyzerHelper
