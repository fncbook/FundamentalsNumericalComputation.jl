"""
horner(c,x)

Evaluate a polynomial whose coefficients are given in ascending
order in `c`, at the point `x`, using Horner's rule.
"""
function horner(c,x)

n = length(c)
y = c[n]
for k in n-1:-1:1
    y = x*y + c[k]
end

return y
end
