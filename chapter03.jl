"""
lsnormal(A,b)

Solve a linear least squares problem by the normal equations.
Returns the minimizer of ||b-Ax||.
"""
function lsnormal(A,b)

N = A'*A;  z = A'*b;
R = cholesky(N).U
w = forwardsub(R',z)                   # solve R'z=c
x = backsub(R,w)                       # solve Rx=z

return x
end


"""
lsqrfact(A,b)

Solve a linear least squares problem by QR factorization. Returns
the minimizer of ||b-Ax||.
"""
function lsqrfact(A,b)

Q,R = qr(A)
c = Q'*b
x = backsub(R,c)

return x
end

"""
qrfact(A)

QR factorization by Householder reflections. Returns Q and R.
"""
function qrfact(A)

    m,n = size(A)
    Qt = Matrix(Diagonal(ones(m)))
    R = float(copy(A))
    for k in 1:n
      z = R[k:m,k]
      v = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
      nrmv = norm(v)
      if nrmv < eps() continue; end  # skip this iteration
      v = v / nrmv;                  # simplifies other formulas
      # Apply the reflection to each relevant column of A and Q
      for j in 1:n
        R[k:m,j] -= v*( 2*(v'*R[k:m,j]) )
      end
      for j in 1:m
        Qt[k:m,j] -= v*( 2*(v'*Qt[k:m,j]) )
      end
    end
    
    return Qt',triu(R)
end
    