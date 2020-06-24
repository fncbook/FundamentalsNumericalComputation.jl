"""
newton(f,dfdx,x1)

Use Newton's method to find a root of `f` starting from `x1`, where
`dfdx` is the derivative of `f`. Returns a vector of root estimates.
"""
function newton(f,dfdx,x1)
    # Operating parameters.
    funtol = 100*eps();  xtol = 100*eps();  maxiter = 40;

    x = [x1]
    y = f(x1)
    dx = Inf   # for initial pass below
    k = 1

    while (abs(dx) > xtol) && (abs(y) > funtol) && (k < maxiter)
        dydx = dfdx(x[k])
        dx = -y/dydx            # Newton step
        push!(x,x[k]+dx)        # append new estimate

        k = k+1
        y = f(x[k])
    end

    if k==maxiter
        @warn "Maximum number of iterations reached."
    end

    return x
end

"""
secant(f,x1,x2)

Use the secant method to find a root of `f` starting from `x1` and
`x2`. Returns a vector of root estimates.
"""
function secant(f,x1,x2)
    # Operating parameters.
    funtol = 100*eps();  xtol = 100*eps();  maxiter = 40;

    x = [x1,x2]
    y1 = f(x1); y2 = 100;
    dx = Inf   # for initial pass below
    k = 2

    while (abs(dx) > xtol) && (abs(y2) > funtol) && (k < maxiter)
        y2 = f(x[k])
        dx = -y2 * (x[k]-x[k-1]) / (y2-y1)   # secant step
        push!(x,x[k]+dx)        # append new estimate

        k = k+1
        y1 = y2    # current f-value becomes the old one next time
    end

    if k==maxiter
        @warn "Maximum number of iterations reached."
    end

    return x
end

"""
newtonsys(f,jac,x1)

Use Newton's method to find a root of a system of equations,
starting from `x1`. The functions `f` and `jac should return the
residual vector and the Jacobian matrix, respectively. Returns
history of root estimates as a vector of vectors.
"""
function newtonsys(f,jac,x1)
    # Operating parameters.
    funtol = 1000*eps();  xtol = 1000*eps();  maxiter = 40;

    x = [float(x1)]
    y,J = f(x1),jac(x1)
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
fdjac(f,x0,y0)

Compute a finite-difference approximation of the Jacobian matrix for
`f` at `x0`, where `y0`=`f(x0)` is given.
"""
function fdjac(f,x0,y0)

delta = sqrt(eps())   # FD step size
m,n = length(y0),length(x0)
if n==1
    J = (f(x0+delta) - y0) / delta
else
    J = zeros(m,n)
    In = I(n)
    for j = 1:n
        J[:,j] = (f(x0 + delta*In[:,j]) - y0) / delta
    end
end

return J
end

"""
levenberg(f,x1,tol)

Use Levenberg's quasi-Newton iteration to find a root of the system
`f`, starting from `x1`, with `tol` as the stopping tolerance in
both step size and residual norm. Returns root estimates as a
matrix, one estimate per column.
"""
function levenberg(f,x1,tol=1e-12)

# Operating parameters.
ftol = tol;  xtol = tol;  maxiter = 40;

x = zeros(length(x1),maxiter)
x = [float(x1)]
fk = f(x1)
k = 1;  s = Inf;
Ak = fdjac(f,x1,fk)   # start with FD Jacobian
jac_is_new = true

lambda = 10;
while (norm(s) > xtol) && (norm(fk) > ftol) && (k < maxiter)
    # Compute the proposed step.
    B = Ak'*Ak + lambda*I
    z = Ak'*fk
    s = -(B\z)

    xnew = x[k] + s
    fnew = f(xnew)

    # Do we accept the result?
    if norm(fnew) < norm(fk)    # accept
        y = fnew - fk
        push!(x,xnew)
        fk = fnew
        k += 1

        lambda = lambda/10   # get closer to Newton
        # Broyden update of the Jacobian.
        Ak = Ak + (y-Ak*s)*(s'/(s'*s))
        jac_is_new = false
    else                       # don't accept
        # Get closer to steepest descent.
        lambda = 4lambda
        # Re-initialize the Jacobian if it's out of date.
        if !jac_is_new
            Ak = fdjac(f,x[k],fk)
            jac_is_new = true
        end
    end
end

if (norm(fk) > 1e-3)
    @warn "Iteration did not find a root."
end

return x
end
