"""
    forwardsub(L,b)

Solve the lower-triangular linear system with matrix `L` and
right-hand side vector `b`.
"""
function forwardsub(L,b)

n = size(L,1)
x = zeros(n)
x[1] = b[1]/L[1,1]
for i in 2:n
    s = sum( L[i,j]*x[j] for j in 1:i-1 )
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
for i in n-1:-1:1
    s = sum( U[i,j]*x[j] for j in i+1:n )
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
L = diagm(ones(n))   # ones on main diagonal, zeros elsewhere
U = zeros(n,n)
Aₖ = float(copy(A))

# Reduction by outer products
for k in 1:n-1
    U[k,:] = Aₖ[k,:]
    L[:,k] = Aₖ[:,k]/U[k,k]
    Aₖ -= L[:,k]*U[k,:]'
end
U[n,n] = Aₖ[n,n]

return L,U
end

"""
    plufact(A)

Compute the PLU factorization of square matrix `A`, returning the
triangular factors and a row permutation vector.
"""
function plufact(A)

n = size(A,1)
L = diagm(ones(n))   # ones on main diagonal, zeros elsewhere
U = zeros(n,n)
p = fill(0,n)
Aₖ = float(copy(A))

# Reduction by outer products
for k in 1:n-1
    p[k] = argmax(abs.(Aₖ[:,k]))
    U[k,:] = Aₖ[p[k],:]
    L[:,k] = Aₖ[:,k]/U[k,k]
    Aₖ -= L[:,k]*U[k,:]'
end
p[n] = argmax(abs.(Aₖ[:,n]))
U[n,n] = Aₖ[p[n],n]

return L[p,:],U,p
end