# Fundamentals of Numerical Computation

These are core functions for the text *Fundamentals of Numerical Computation* by T. A. Driscoll and R. J. Braun ([preview site](http://tobydriscoll.net/unlinked/fnc-preview/)). They are a companion to the Julia (2nd) edition, to be published in 2022.

**For Julia versions of the functions accompanying the MATLAB (first) edition, go to https://github.com/fncbook/fnc instead.**

## Installation

1. [Install Julia.](https://julialang.org/downloads/) Any of the available methods for the latest stable version 1.x should be fine. The Julia Pro version comes with an integrated editor and many preinstalled packages, and it might be the best choice for those not comfortable with their command line shell.
2. Start Julia on your machine. Look for the `julia>` prompt.
3. At the prompt type the `]` character. This will turn the prompt a different color and say `pkg>`, indicating that you can give commands to the package manager.
4. At the `pkg` prompt, type

   ```julia
   add https://github.com/fncbook/FundamentalsNumericalComputation.jl
   ```

   The process will take a while. In order to flatten the learning curve, this package loads many other standard numerical and plotting packages and makes them available, so there is a lot of code to install and compile.
5. Hit the backspace (delete on Mac) key to go back to the main Julia prompt. **Steps 3-5 should only need to be done once per Julia installation.**
6. In order to use the functions, in each new Julia session you must enter

   ```julia
   using FundamentalsNumericalComputation
   ```

   After this completes, all the text's functions can be accessed with the prefix `FNC`. E.g., `FNC.lufact`, `FNC.rk23`, etc.

## Usage

None of the functions are exported, so they must all be prefixed with the package name. However, the constant `FNC` is set as an alias to the package name to make typing more reasonable. 

There is a bare-bones [documentation site](https://fncbook.github.io/FundamentalsNumericalComputation.jl/functions/), but it only gives summaries of the documentation strings. The textbook is meant to be the real guide.

### Startup speed

After installation, importing the package with a `using` statement should only take 5-10 seconds in Julia 1.6 or later. The first plot created in a Julia session may take 20-60 seconds to appear. This is a well-known irritant in the Julia ecosystem, and progress is steadily being made toward reducing the lag. In the meantime, power users might investigate using [`PackageCompiler`](https://julialang.github.io/PackageCompiler.jl/dev/) to get better startup performance.

## Alternatives 

These codes are for instructional purposes. They are not recommended for applications. Compared to superior alternatives, they lack generality, efficiency, and robustness to syntax and more subtle mistakes. My personal recommendations for preferable alternatives are as follows. Most of them are demonstrated in the textbook.

|  Problem type | Associated Julia packages |
|-----------------|---------------|
| Linear system / least squares   |  [`LinearAlgebra`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#man-linalg) (standard library)   |
| Sparse matrix     | [SparseArrays](https://docs.julialang.org/en/v1/stdlib/SparseArrays/#Sparse-Arrays), [IterativeSolvers](https://iterativesolvers.julialinearalgebra.org/stable/),  [Arpack](https://arpack.julialinearalgebra.org/stable/), [Preconditioners](https://github.com/mohamed82008/Preconditioners.jl) |
| Polynomial interpolation/approximation  | [Polynomials](https://juliamath.github.io/Polynomials.jl/stable/), [ApproxFun](https://juliaapproximation.github.io/ApproxFun.jl/stable/) |
| Polynomial roots    | [Polynomials](https://juliamath.github.io/Polynomials.jl/stable/#Root-finding-1) |
| Rootfinding       | [NLsolve](https://github.com/JuliaNLSolvers/NLsolve.jl) |
| Finite differences | [FiniteDifferences](https://juliadiff.org/FiniteDifferences.jl/latest/), [FiniteDiff](https://github.com/JuliaDiff/FiniteDiff.jl) | 
| Integration        | [QuadGK](https://juliamath.github.io/QuadGK.jl/stable/) |
| Spline  | [Dierckx](https://github.com/kbarbary/Dierckx.jl)  | 
| [Initial-value problem](https://diffeq.sciml.ai/latest/tutorials/ode_example/#ode_example)   | [DifferentialEquations](https://diffeq.sciml.ai/latest/) | 
| [Boundary-value problem](https://diffeq.sciml.ai/latest/tutorials/bvp_example/#Boundary-Value-Problems)  | [DifferentialEquations](https://diffeq.sciml.ai/latest/) |
| Method of lines  |  [MethodOfLines](https://methodoflines.sciml.ai/dev/) (still under development)



## License

This code generating the site is under an MIT license. Please see the LICENSE file for details.
