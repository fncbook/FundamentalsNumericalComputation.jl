"""
    ndgrid(x,y,...)

Given ``d`` vector inputs, returns ``d`` matrices representing the coordinate
functions on the tensor product grid.
"""
function ndgrid(x...)
    I = CartesianIndices( fill(undef,length.(x)) )
    return [ [ x[d][i[d]] for i in I]  for d in 1:length(x) ]
end

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
    x,Dx,Dxx = FNC.diffcheb(m,xspan)
    y,Dy,Dyy = FNC.diffcheb(n,yspan)
    mtx = h -> [ h(x,y) for x in x, y in y ]
    X,Y = mtx((x,y)->x),mtx((x,y)->y)
    N = (m+1)*(n+1)   # total number of unknowns

    # Form the collocated PDE as a linear system.
    A = kron(I(n+1),Dxx) + kron(Dyy,I(m+1))  # Laplacian matrix
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
    elliptic(f,g,m,xspan,n,yspan)

Solve the elliptic PDE `f`(x,y,u,u_x,u_xx,u_y,u_yy)=0 on the 
rectangle `xspan` x `yspan`, subject to `g`(x,y)=0 on the boundary. 
Uses `m`+1 points in x by `n`+1 points in y in a Chebyshev 
discretization.

Returns vectors defining the grid and a matrix of grid solution values.
"""
function elliptic(f,g,m,xspan,n,yspan)
    # Discretize the domain.
    x,Dx,Dxx = FNC.diffcheb(m,xspan)
    y,Dy,Dyy = FNC.diffcheb(n,yspan)
    mtx = h -> [ h(x,y) for x in x, y in y ]
    unvec = u -> reshape(u,m+1,n+1)
    N = (m+1)*(n+1)   # total number of unknowns
    X,Y = mtx((x,y)->x),mtx((x,y)->y)
    
    # Identify boundary locations and evaluate the boundary condition.
    isboundary = trues(m+1,n+1)
    isboundary[2:m,2:n] .= false
    idx = vec(isboundary)
    gb = g.(X[idx],Y[idx])

    # Evaluate the discretized PDE and its Jacobian, with all the
    # boundary condition modifications applied.
    function residual(u)
        U = unvec(u)
        R = f(x,y,U,Dx*U,Dxx*U,U*Dy',U*Dyy')
        @. R[idx] = u[idx] - gb
        return vec(R)
    end

    # Solve the equation.
    u = levenberg(residual,zeros(N))[end]
    U = unvec(u)

    return function (ξ,η)
		v = [ chebinterp(x,u,ξ) for u in eachcol(U) ]
		return chebinterp(y,v,η)
    end

end

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
