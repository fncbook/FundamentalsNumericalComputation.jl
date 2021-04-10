"""
    newton(f,dfdx,x₁)

Use Newton's method to find a root of `f` starting from `x₁`, where
`dfdx` is the derivative of `f`. Returns a vector of root estimates.
"""
function newton(f,dfdx,x₁)

# Operating parameters.
funtol = 100*eps();  xtol = 100*eps();  maxiter = 40;

x = [x₁]
y = f(x₁)
dx = Inf   # for initial pass below
k = 1

while (abs(dx) > xtol) && (abs(y) > funtol) && (k < maxiter)
    dydx = dfdx(x[k])
    dx = -y/dydx            # Newton step
    push!(x,x[k]+dx)        # append new estimate

    k += 1
    y = f(x[k])
end

if k==maxiter
    @warn "Maximum number of iterations reached."
end

return x
end

"""
    secant(f,x₁,x₂)

Use the secant method to find a root of `f` starting from `x₁` and
`x₂`. Returns a vector of root estimates.
"""
function secant(f,x₁,x₂)

# Operating parameters.
funtol = 100*eps();  xtol = 100*eps();  maxiter = 40;

x = [x₁,x₂]
y₁ = f(x₁); y₂ = 100;
dx = Inf   # for initial pass below
k = 2

while (abs(dx) > xtol) && (abs(y₂) > funtol) && (k < maxiter)
    y₂ = f(x[k])
    dx = -y₂ * (x[k]-x[k-1]) / (y₂-y₁)   # secant step
    push!(x,x[k]+dx)        # append new estimate

    k += 1
    y₁ = y₂    # current f-value becomes the old one next time
end

if k==maxiter
    @warn "Maximum number of iterations reached."
end

return x
end

"""
    newtonsys(f,jac,x₁)

Use Newton's method to find a root of a system of equations,
starting from `x₁`. The functions `f` and `jac should return the
residual vector and the Jacobian matrix, respectively. Returns the
history of root estimates as a vector of vectors.
"""
function newtonsys(f,jac,x₁)

# Operating parameters.
funtol = 1000*eps();  xtol = 1000*eps();  maxiter = 40;

x = [float(x₁)]
y,J = f(x₁),jac(x₁)
dx = Inf   # for initial pass below
k = 1

while (norm(dx) > xtol) && (norm(y) > funtol) && (k < maxiter)
    dx = -(J\y)             # Newton step
    push!(x,x[k] + dx)    # append to history
    k += 1
    y,J = f(x[k]),jac(x[k])
end

if k==maxiter
    @warn "Maximum number of iterations reached."
end

return x
end

"""
    fdjac(f,x₀,y₀)

Compute a finite-difference approximation of the Jacobian matrix for
`f` at `x₀`, where `y₀`=`f(x₀)` is given.
"""
function fdjac(f,x₀,y₀)

δ = sqrt(eps())   # FD step size
m,n = length(y₀),length(x₀)
if n==1
    J = (f(x₀+δ) - y₀) / δ
else
    J = zeros(m,n)
    for j in 1:n
        x = copy(x₀)
        x[j] += δ
        J[:,j] = (f(x) - y₀) / δ
    end
end

return J
end

"""
    levenberg(f,x₁,tol)

Use Levenberg's quasi-Newton iteration to find a root of the system
`f`, starting from `x₁`, with `tol` as the stopping tolerance in
both step size and residual norm. Returns the history of root estimates 
as a vector of vectors.
"""
function levenberg(f,x₁,tol=1e-12)

# Operating parameters.
ftol = tol;  xtol = tol;  maxiter = 40;

x = zeros(length(x₁),maxiter)
x = [float(x₁)]
fₖ = f(x₁)
k = 1;  s = Inf;
Aₖ = fdjac(f,x₁,fₖ)   # start with FD Jacobian
jac_is_new = true

λ = 10;
while (norm(s) > xtol) && (norm(fₖ) > ftol) && (k < maxiter)
    # Compute the proposed step.
    B = Aₖ'*Aₖ + λ*I
    z = Aₖ'*fₖ
    s = -(B\z)

    xnew = x[k] + s
    fnew = f(xnew)

    # Do we accept the result?
    if norm(fnew) < norm(fₖ)    # accept
        y = fnew - fₖ
        push!(x,xnew)
        fₖ = fnew
        k += 1

        λ = λ/10   # get closer to Newton
        # Broyden update of the Jacobian.
        Aₖ = Aₖ + (y-Aₖ*s)*(s'/(s'*s))
        jac_is_new = false
    else                       # don't accept
        # Get closer to steepest descent.
        λ = 4λ
        # Re-initialize the Jacobian if it's out of date.
        if !jac_is_new
            Aₖ = fdjac(f,x[k],fₖ)
            jac_is_new = true
        end
    end
end

if (norm(fₖ) > 1e-3)
    @warn "Iteration did not find a root."
end

return x
end
