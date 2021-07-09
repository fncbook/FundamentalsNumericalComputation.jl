"""
    diffper(n,xspan)

Construct 2nd-order differentiation matrices for functions with
periodic end conditions, using `n` unique nodes in the interval
`xspan`. Returns a vector of nodes and the matrices for the first
and second derivatives.
"""
function diffper(n,xspan)
    a,b = xspan
    h = (b-a)/n
    x = @. a + h*(0:n-1)   # nodes, omitting the repeated data

    # Construct Dx by diagonals, then correct the corners.
    dp = fill(0.5/h,n-1)        # superdiagonal
    dm = fill(-0.5/h,n-1)       # subdiagonal
    Dx = diagm(-1=>dm,1=>dp)
    Dx[1,n] = -1/(2*h)
    Dx[n,1] = 1/(2*h)

    # Construct Dxx by diagonals, then correct the corners.
    d0 =  fill(-2/h^2,n)        # main diagonal
    dp =  ones(n-1)/h^2         # superdiagonal and subdiagonal
    Dxx = diagm(-1=>dp,0=>d0,1=>dp)
    Dxx[1,n] = 1/(h^2)
    Dxx[n,1] = 1/(h^2)

    return x,Dx,Dxx
end

"""
    parabolic(ϕ,xspan,m,g₁,g₂,tspan,init)

Solve a parabolic PDE by the method of lines. The PDE is 
∂u/∂t = `ϕ`(t,x,u,∂u/∂x,∂^2u/∂x^2), `xspan` gives the space 
domain, m gives the degree of a Chebyshev spectral discretization, 
`g₁` and `g₂` are functions of (u,∂u/∂x) at the domain ends that 
should be made zero, `tspan` is the time domain, and `init` is a 
function of x that gives the initial condition. Returns a vector
`x` and a function of t that gives the semidiscrete solution at `x`. 
"""
function parabolic(ϕ,xspan,m,g₁,g₂,tspan,init)
    x,Dₓ,Dₓₓ = diffcheb(m,xspan)
    int = 2:m    # indexes of interior nodes

    function extend(v)
        function objective(ubc)
            u₀,uₘ = ubc
            uₓ = Dₓ*[u₀;v;uₘ]
            return [g₁(u₀,uₓ[1]),g₂(uₘ,uₓ[end])]
        end
        ubc = levenberg(objective,[0,0])[end]
        return [ubc[1];v;ubc[2]]
    end

    function ode!(f,v,p,t)
        u = extend(v)
        uₓ,uₓₓ = Dₓ*u,Dₓₓ*u
        f .= ϕ(t,x[int],u[int],uₓ[int],uₓₓ[int])
    end

    ivp = ODEProblem(ode!,init.(x[int]),float.(tspan))
    u = solve(ivp)

    return x,t->extend(u(t))
end