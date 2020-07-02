"""
forwardsub(L,b)

Solve the lower-triangular linear system with matrix `L` and
right-hand side vector `b`.
"""
function forwardsub(L,b)

n = size(L,1)
x = zeros(n)
x[1] = b[1]/L[1,1]
for i = 2:n
    s = sum( L[i,j]*x[j] for j=1:i-1 )
    x[i] = ( b[i] - s ) / L[i,i]
end

return x
end

"""
backsub(U,b)

Solve the upper-triangular linear system with matrix `U` and
right-hand side vector `b`.
"""
function backsub(U,b)

n = size(U,1)
x = zeros(n)
x[n] = b[n]/U[n,n]
for i = n-1:-1:1
    s = sum( U[i,j]*x[j] for j=i+1:n )
    x[i] = ( b[i] - s ) / U[i,i]
end

return x
end

"""
lufact(A)

Compute the LU factorization of square matrix `A`, returning the
factors.
"""
function lufact(A)

n = size(A,1)
L = diagm(0=>ones(n))  # ones on main diagonal, zeros elsewhere
U = float(copy(A))

# Gaussian elimination
for j = 1:n-1
  for i = j+1:n
    L[i,j] = U[i,j] / U[j,j]   # row multiplier
    U[i,j:n] -= L[i,j]*U[j,j:n]
  end
end

return L,triu(U)
end
