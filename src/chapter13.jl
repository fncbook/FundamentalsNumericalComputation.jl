"""
    poisson(f,g,m,xspan,n,yspan)

Solve Poisson's equation on a rectangle by finite differences.
Function `f` is the forcing function and function `g` gives the
Dirichlet boundary condition. The rectangle is the tensor product of
intervals `xspan` and `yspan`,  and the discretization uses `m`+1
and `n`+1 points in the two coordinates.

Returns vectors defining the grid and a matrix of grid solution values.
"""
function poisson(f,g,m,xspan,n,yspan)
    # Discretize the domain.
    x,Dx,Dxx = FNC.diffmat2(m,xspan)
    y,Dy,Dyy = FNC.diffmat2(n,yspan)
    mtx = h -> [ h(x,y) for x in x, y in y ]
    X,Y = mtx((x,y)->x),mtx((x,y)->y)
    N = (m+1)*(n+1)   # total number of unknowns

    # Form the collocated PDE as a linear system.
    A = kron(I(n+1),sparse(Dxx)) + kron(sparse(Dyy),I(m+1))
    b = vec( mtx(f) )

    # Identify boundary locations.
    isboundary = trues(m+1,n+1)
    isboundary[2:m,2:n] .= false
    idx = vec(isboundary)

    # Apply Dirichlet condition.
    scale = maximum(abs,A[n+2,:])
    A[idx,:] = scale * I(N)[idx,:]        # Dirichet assignment
    b[idx] = scale * g.(X[idx],Y[idx])    # assigned values

    # Solve the linear sytem and reshape the output.
    u = A\b
    U = reshape(u,m+1,n+1)
    return x,y,U
end

"""
    elliptic(ϕ,g,m,xspan,n,yspan)

Solve the elliptic PDE `ϕ`(x,y,u,u_x,u_xx,u_y,u_yy)=0 on the 
rectangle `xspan` x `yspan`, subject to `g`(x,y)=0 on the boundary. 
Uses `m`+1 points in x by `n`+1 points in y in a Chebyshev 
discretization.

Returns vectors defining the grid and a matrix of grid solution values.
"""
function elliptic(ϕ,g,m,xspan,n,yspan)
    # Discretize the domain.
    x,Dx,Dxx = diffcheb(m,xspan)
    y,Dy,Dyy = diffcheb(n,yspan)
    X = [ x for x in x, y in y ]
    Y = [ y for x in x, y in y ]
    unvec = u -> reshape(u,m+1,n+1)
    
    # Identify boundary locations and evaluate the boundary condition.
    isboundary = trues(m+1,n+1)
    isboundary[2:m,2:n] .= false
    idx = vec(isboundary)
    gb = g.(X[idx],Y[idx])

    # Evaluate the discretized PDE and its Jacobian, with all the
    # boundary condition modifications applied.
    function residual(u)
        U = unvec(u)
        R = ϕ(X,Y,U,Dx*U,Dxx*U,U*Dy',U*Dyy')
        @. R[idx] = u[idx] - gb
        return vec(R)
    end

    # Solve the equation.
    u = levenberg(residual,vec(zeros(size(X))))[end]
    U = unvec(u)

    return function (ξ,η)
		v = [ chebinterp(x,u,ξ) for u in eachcol(U) ]
		return chebinterp(y,v,η)
    end

end

"Evaluate Chebyshev interpolant with nodes x, values v, at point ξ"
function chebinterp(x,v,ξ)
    n = length(x)-1
    w = (-1.0).^(0:n)
    w[[1,n+1]] .*= 0.5

	terms = @. w / (ξ - x)
	if any(isinf.(terms))     # exactly at a node
		idx = findfirst(ξ.==x)
		return v[idx]
	else
		return sum(v.*terms) / sum(terms)
	end
end
