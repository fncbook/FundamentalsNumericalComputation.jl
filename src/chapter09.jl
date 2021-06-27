"""
    polyinterp(t,y)

Construct a callable polynomial interpolant through the points in
vectors `t`,`y` using the barycentric interpolation formula.
"""
function polyinterp(t,y)
    n = length(t)-1
    C = (t[n+1]-t[1]) / 4           # scaling factor to ensure stability
    tc = t/C

    # Adding one node at a time, compute inverses of the weights.
    ω = ones(n+1)
    for m in 0:n-1
        d = tc[1:m+1] .- tc[m+2]    # vector of node differences
        @. ω[1:m+1] *= d            # update previous
        ω[m+2] = prod( -d )         # compute the new one
    end
    w = 1 ./ ω                      # go from inverses to weights

    # This function evaluates the interpolant at given x.
    p = function (x)
        terms = @. w / (x - t)
        if any(isinf.(terms))     # there was division by zero
            # return the node's data value
            idx = findfirst(x.==t)
            f = y[idx]
        else
            f = sum(y.*terms) / sum(terms)
        end
    end
    return p
end

"""
    triginterp(t,y)

Construct the trigonometric interpolant for the points defined by 
vectors `t` and `y`.
"""
function triginterp(t,y)
    N = length(t)

    function τ(x)
        if x==0
            return 1.0
        else
            denom = isodd(N) ? N*sin(π*x/2) : N*tan(π*x/2)
            return sin(N*π*x/2)/denom
        end
    end

    return function (x)
        sum( y[k]*τ(x-t[k]) for k in eachindex(y) )
    end
end

"""
    ccint(f,n)

Perform Clenshaw-Curtis integration for the function `f` on `n`+1
nodes in [-1,1]. Returns the integral estimate and a vector of the 
nodes used. Note: `n` must be even.
"""
function ccint(f,n)
    @assert iseven(n) "Value of `n` must be an even integer."
    # Find Chebyshev extreme nodes.
    θ = [ i*π/n for i in 0:n ]
    x = -cos.(θ)

    # Compute the C-C weights.
    c = similar(θ)
    c[[1,n+1]] .= 1/(n^2-1)
    s = sum( cos.(2k*θ[2:n])/(4k^2-1) for k in 1:n/2-1 )
    v = @. 1 - 2s - cos(n*θ[2:n])/(n^2-1)
    c[2:n] = 2v/n

    # Evaluate integrand and integral.
    I = dot(c,f.(x))   # vector inner product
    return I,x
end

"""
    glint(f,n)

Perform Gauss-Legendre integration for the function `f` on `n` nodes
in (-1,1). Returns the integral estimate and a vector of the nodes used.
"""
function glint(f,n)
    # Nodes and weights are found via a tridiagonal eigenvalue problem.
    β = @. 0.5/sqrt(1-(2*(1:n-1))^(-2))
    T = diagm(-1=>β,1=>β)
    λ,V = eigen(T)
    p = sortperm(λ)
    x = λ[p]               # nodes
    c = @. 2V[1,p]^2       # weights

    # Evaluate the integrand and compute the integral.
    I = dot(c,f.(x))      # vector inner product
    return I,x
end

"""
    intinf(f,tol)

Perform adaptive doubly-exponential integration of function `f` 
over (-Inf,Inf), with error tolerance `tol`. Returns the integral 
estimate and a vector of the nodes used.
"""
function intinf(f,tol)   
    x = t -> sinh(sinh(t))
    dx_dt = t -> cosh(t)*cosh(sinh(t))
    g = t -> f(x(t))*dx_dt(t)

    # Find where to truncate the integration interval.
    M = 3
    while (abs(g(-M)) > tol/100) || (abs(g(M)) > tol/100)
        M += 0.5
        if M > log(log(floatmax()))
            @warn "Function may not decay fast enough."
            M -= 0.5
            break
        end
    end

    I,t = intadapt(g,-M,M,tol)
	return I,x.(t)
end

"""
    intsing(f,tol)

Adaptively integrate function `f` over (0,1), where `f` may be 
singular at zero, with error tolerance `tol`. Returns the
integral estimate and a vector of the nodes used.
"""
function intsing(f,tol)
    x = t -> 2/(1+exp(2sinh(t)))
	dx_dt = t -> cosh(t)/cosh(sinh(t))^2
	g = t -> f(x(t))*dx_dt(t)

    # Find where to truncate the integration interval.
    M = 3
    while abs(g(M)) > tol/100
        M += 0.5
        if M > log(log(2/floatmin()))
            @warn "Function may grow too rapidly."
            M -= 0.5
            break
        end
    end

    I,t = intadapt(g,0,M,tol)
	return I,x.(t)
end
