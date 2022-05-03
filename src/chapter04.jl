"""
    newton(f,dfdx,x₁[;maxiter,ftol,xtol])

Use Newton's method to find a root of `f` starting from `x₁`, where
`dfdx` is the derivative of `f`. Returns a vector of root estimates.

The optional keyword parameters set the maximum number of iterations
and the stopping tolerance for values of `f` and changes in `x`.
"""
function newton(f,dfdx,x₁;maxiter=40,ftol=100*eps(),xtol=100*eps())
    x = [float(x₁)]
    y = f(x₁)
    Δx = Inf   # for initial pass below
    k = 1

    while (abs(Δx) > xtol) && (abs(y) > ftol)
        dydx = dfdx(x[k])
        Δx = -y/dydx            # Newton step
        push!(x,x[k]+Δx)        # append new estimate

        k += 1
        y = f(x[k])
        if k==maxiter
            @warn "Maximum number of iterations reached."
            break   # exit loop
        end
    end
    return x
end

"""
    secant(f,x₁,x₂[;maxiter,ftol,xtol])

Use the secant method to find a root of `f` starting from `x₁` and
`x₂`. Returns a vector of root estimates.

The optional keyword parameters set the maximum number of iterations
and the stopping tolerance for values of `f` and changes in `x`.
"""
function secant(f,x₁,x₂;maxiter=40,ftol=100*eps(),xtol=100*eps())
    x = [float(x₁),float(x₂)]
    y₁ = f(x₁)
    Δx,y₂ = Inf,Inf   # for initial pass in the loop below
    k = 2

    while (abs(Δx) > xtol) && (abs(y₂) > ftol) 
        y₂ = f(x[k])
        Δx = -y₂ * (x[k]-x[k-1]) / (y₂-y₁)   # secant step
        push!(x,x[k]+Δx)        # append new estimate

        k += 1
        y₁ = y₂    # current f-value becomes the old one next time
        
        if k==maxiter
            @warn "Maximum number of iterations reached."
            break   # exit loop
        end
    end
    return x
end

"""
    newtonsys(f,jac,x₁[;maxiter,ftol,xtol])

Use Newton's method to find a root of a system of equations,
starting from `x₁`. The functions `f` and `jac` should return the
residual vector and the Jacobian matrix, respectively. Returns the
history of root estimates as a vector of vectors.

The optional keyword parameters set the maximum number of iterations
and the stopping tolerance for values of `f` and changes in `x`.

"""
function newtonsys(f,jac,x₁;maxiter=40,ftol=1000*eps(),xtol=1000*eps())
    x = [float(x₁)]
    y,J = f(x₁),jac(x₁)
    Δx = Inf   # for initial pass below
    k = 1

    while (norm(Δx) > xtol) && (norm(y) > ftol)
        Δx = -(J\y)             # Newton step
        push!(x,x[k] + Δx)    # append to history
        k += 1
        y,J = f(x[k]),jac(x[k])

        if k==maxiter
            @warn "Maximum number of iterations reached."
            break
        end
    end
    return x
end

"""
    fdjac(f,x₀[,y₀])

Compute a finite-difference approximation of the Jacobian matrix for
`f` at `x₀`, where `y₀`=`f(x₀)` may be given.
"""
function fdjac(f,x₀,y₀=f(x₀))
    δ = sqrt(eps())*max(norm(x₀),1)   # FD step size
    m,n = length(y₀),length(x₀)
    if n==1
        J = (f(x₀+δ) - y₀) / δ
    else
        J = zeros(m,n)
        x = copy(x₀)
        for j in 1:n
            x[j] += δ
            J[:,j] = (f(x) - y₀) / δ
            x[j] -= δ
        end
    end
    return J
end

"""
    levenberg(f,x₁[;maxiter,ftol,xtol])

Use Levenberg's quasi-Newton iteration to find a root of the system
`f` starting from `x₁` Returns the history of root estimates as a 
vector of vectors.

The optional keyword parameters set the maximum number of iterations
and the stopping tolerance for values of `f` and changes in `x`.

"""
function levenberg(f,x₁;maxiter=40,ftol=1e-12,xtol=1e-12)
    x = [float(x₁)]
    yₖ = f(x₁)
    k = 1;  s = Inf;
    A = fdjac(f,x[k],yₖ)   # start with FD Jacobian
    jac_is_new = true

    λ = 10;
    while (norm(s) > xtol) && (norm(yₖ) > ftol)
        # Compute the proposed step.
        B = A'*A + λ*I
        z = A'*yₖ
        s = -(B\z)
        
        x̂ = x[k] + s
        ŷ = f(x̂)

        # Do we accept the result?
        if norm(ŷ) < norm(yₖ)    # accept
            λ = λ/10   # get closer to Newton
            # Broyden update of the Jacobian.
            A += (ŷ-yₖ-A*s)*(s'/(s'*s))
            jac_is_new = false
            
            push!(x,x̂)
            yₖ = ŷ
            k += 1
        else                       # don't accept
            # Get closer to gradient descent.
            λ = 4λ
            # Re-initialize the Jacobian if it's out of date.
            if !jac_is_new
                A = fdjac(f,x[k],yₖ)
                jac_is_new = true
            end
        end

        if k==maxiter
            @warn "Maximum number of iterations reached."
            break
        end
        
    end
    return x
end
