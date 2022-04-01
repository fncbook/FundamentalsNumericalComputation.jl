"""
    poweriter(A,numiter)

Perform `numiter` power iterations with the matrix `A`, starting
from a random vector. Returns a vector of eigenvalue estimates
and the final eigenvector approximation.
"""
function poweriter(A,numiter)
    n = size(A,1)
    x = normalize(randn(n),Inf)
    β = zeros(numiter)
    for k in 1:numiter
        y = A*x
        m = argmax(abs.(y))
        β[k] = y[m]/x[m]
        x = y/y[m]
    end
    return β,x
end

"""
    inviter(A,s,numiter)

Perform `numiter` inverse iterations with the matrix `A` and shift
`s`, starting from a random vector. Returns a vector of
eigenvalue estimates and the final eigenvector approximation.
"""
function inviter(A,s,numiter)
    n = size(A,1)
    x = normalize(randn(n),Inf)
    β = zeros(numiter)
    fact = lu(A - s*I)
    for k in 1:numiter
        y = fact\x
        normy,m = findmax(abs.(y))
        β[k] = x[m]/y[m] + s
        x = y/y[m]
    end
    return β,x
end

"""
    arnoldi(A,u,m)

Perform the Arnoldi iteration for `A` starting with vector `u`, out
to the Krylov subspace of degree `m`. Returns the orthonormal basis
(`m`+1 columns) and the upper Hessenberg `H` of size `m`+1 by `m`.
"""
function arnoldi(A,u,m)
    n = length(u)
    Q = zeros(n,m+1)
    H = zeros(m+1,m)
    Q[:,1] = u/norm(u)
    for j in 1:m
        # Find the new direction that extends the Krylov subspace.
        v = A*Q[:,j]
        # Remove the projections onto the previous vectors.
        for i in 1:j
            H[i,j] = dot(Q[:,i],v)
            v -= H[i,j]*Q[:,i]
        end
        # Normalize and store the new basis vector.
        H[j+1,j] = norm(v)
        Q[:,j+1] = v/H[j+1,j]
    end
    return Q,H
end

"""
    gmres(A,b,m)

Do `m` iterations of GMRES for the linear system `A`*x=`b`. Returns
the final solution estimate x and a vector with the history of
residual norms. (This function is for demo only, not practical use.)
"""
function gmres(A,b,m)
    n = length(b)
    Q = zeros(n,m+1)
    Q[:,1] = b/norm(b)
    H = zeros(m+1,m)

    # Initial solution is zero.
    x = 0
    residual = [norm(b);zeros(m)]

    for j in 1:m
        # Next step of Arnoldi iteration.
        v = A*Q[:,j]
        for i in 1:j
            H[i,j] = dot(Q[:,i],v)
            v -= H[i,j]*Q[:,i]
        end
        H[j+1,j] = norm(v)
        Q[:,j+1] = v/H[j+1,j]

        # Solve the minimum residual problem.
        r = [norm(b); zeros(j)]
        z = H[1:j+1,1:j] \ r
        x = Q[:,1:j]*z
        residual[j+1] = norm( A*x - b )
    end
    return x,residual
end

### These are not part of the text proper, but help to construct 
### interesting examples.

"""
    sprandsym(n,density,λ)
    sprandsym(n,density,rcond)

Construct a randomized `n` by `n` symmetric sparse matrix of approximate 
density `density`. For vector `λ`, the matrix has eigenvalues as 
prescribed by λ. For scalar `rcond`, the matrix has condition number 
equal to 1/`rcond`.
"""
function sprandsym(n,density,rcond::Number)
    λ = [ rcond^(i/(n-1)) for i in 0:n-1 ]
    sprandsym(n,density,λ)
end

function sprandsym(n,density,λ::Vector)
    # Apply a random Jacobi rotation similarity transformation.
    function randjr(A)
        theta = 2π*rand()
        c,s = cos(theta),sin(theta)
        i = j = 0
        while i==j
          i,j = rand(1:n,2)
        end
        A[[i,j],:] = [c s;-s c]*A[[i,j],:]
        A[:,[i,j]] = A[:,[i,j]]*[c -s;s c]
        return A
    end

    targetnz = ceil(min(0.98,density)*n^2)
    A = spdiagm(λ)

    while nnz(A) < targetnz
        A = randjr(A)
    end
    return A
end

"""
    poisson(n)

Construct the finite-difference Laplacian matrix for an `n` by `n` 
interior lattice.
"""
function poisson(n)
    D = spdiagm(-1=>fill(-1,n-1),0=>fill(2,n),1=>fill(-1,n-1)) * (n+1)^2/pi^2
    return kron(D,I(n)) + kron(I(n),D)
end
