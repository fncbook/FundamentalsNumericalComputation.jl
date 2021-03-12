"""
polyinterp(t,y)

Return a callable polynomial interpolant through the points in
vectors `t`,`y`. Uses the barycentric interpolation formula.
"""
function polyinterp(t,y)
    @assert (isa(t,OffsetArray) && isa(y,OffsetArray)) "Vectors must be indexed 0:n"
    n = length(t)-1
    C = (t[n]-t[0]) / 4           # scaling factor to ensure stability
    tc = t/C
    
    # Adding one node at a time, compute inverses of the weights.
    ω = OffsetArray(ones(n+1),0:n)
    for m in 0:n-1
        d = tc[0:m] .- tc[m+1]    # vector of node differences
        @. ω[0:m] *= d            # update previous
        ω[m+1] = prod( -d )       # compute the new one
    end
    w = 1 ./ ω                    # go from inverses to weights

    p = function (x)
        # Compute interpolant.
        terms = @. w / (x - t)
        if any(isinf.(terms))     # there was division by zero
            # Apply L'Hôpital's Rule exactly.
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

Return trigonometric interpolant for points defined by vectors `t`
and `y`.
"""
function triginterp(t,y)
    N = length(t)

    function trigcardinal(x)
        if isodd(N)      # odd
            tau = sin(N*π*x/2) / (N*sin(π*x/2))
        else             # even
            tau = sin(N*π*x/2) / (N*tan(π*x/2))
        end
        if isnan(tau)
            tau = 1
        end
        return tau
    end

    p = function (x)
        sum( y[k]*trigcardinal(x-t[k]) for k in eachindex(y) )
    end
    return p
end

"""
ccint(f,n)

Perform Clenshaw-Curtis integration for the function `f` on `n`+1
nodes in [-1,1]. Return integral and a vector of the nodes used.
Note: `n` must be even.
"""
function ccint(f,n)
    # Find Chebyshev extreme nodes.
    θ = OffsetArray([ i*π/n for i in 0:n ],0:n)
    x = -cos.(θ)

    # Compute the C-C weights.
    c = similar(θ)
    c[[0,n]] .= 1/(n^2-1)
    s = sum( cos.(2k*θ[1:n-1])/(4k^2-1) for k in 1:n/2-1 )
    v = @. 1 - 2s - cos(n*θ[1:n-1])/(n^2-1)
    c[1:n-1] = 2v/n

    # Evaluate integrand and integral.
    I = dot(c,f.(x))   # vector inner product
    return I,x
end

"""
glint(f,n)

Perform Gauss-Legendre integration for the function `f` on `n` nodes
in (-1,1). Return integral and a vector of the nodes used.
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
intde(f,h,M)

Perform doubly-exponential integration of function `f` over
(-Inf,Inf), using discretization size `h` and truncation point `M`.
Return integral and a vector of the nodes used.
"""
function intde(f,h,M)
    # Find where to truncate the trapezoid sum.
    K = ceil( log(4/π*log(2M))/h )

    # Integrate by trapezoids in a transformed variable t.
    t = h*(-K:K)
    x = @. sinh(π/2*sinh(t))
    dxdt = @. π/2*cosh(t)*cosh(π/2*sinh(t))

    I = h*dot(f.(x),dxdt)
    return I,x
end

"""
intsing(f,h,δ)

Integrate function `f` (possibly singular at 1 and -1) over
[-1+`δ`,1-`δ`] using discretization size `h`. Return
integral and a vector of the nodes used.
"""
function intsing(f,h,δ)
    # Find where to truncate the trapezoid sum.
    K = ceil(log(-2/π*log(δ/2))/h)

    # Integrate over a transformed variable.
    t = h*(-K:K)
    x = @. tanh(π/2*sinh(t))
    dxdt = @. π/2*cosh(t) / (cosh(π/2*sinh(t))^2)

    I = h*dot(f.(x),dxdt)
    return I,x
end
