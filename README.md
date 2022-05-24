# Fundamentals of Numerical Computation

[![][docs-stable-img]][docs-stable-url]

These are core functions for the text *Fundamentals of Numerical Computation* by T. A. Driscoll and R. J. Braun ([preview site](http://tobydriscoll.net/unlinked/fnc-preview/)). They are a companion to the Julia (2nd) edition, to be published in 2022.

**For Julia versions of the functions accompanying the MATLAB (first) edition, go to https://github.com/fncbook/fnc instead.**

## Installation

1. [Install Julia.](https://julialang.org/downloads/) Any of the available methods for the latest stable version 1.x should be fine. The Julia Pro version comes with an integrated editor and many preinstalled packages, and it might be the best choice for those not comfortable with their command line shell.
2. Start Julia on your machine. Look for the `julia>` prompt.
3. At the prompt type the `]` character. This will turn the prompt a different color and say `pkg>`, indicating that you can give commands to the package manager.
4. At the `pkg` prompt, type

   ```julia
   add FundamentalsNumericalComputation
   ```

   The process will take a while. In order to flatten the learning curve, this package loads many other standard numerical and plotting packages and makes them available, so there is a lot of code to install and compile.
5. Hit the backspace (delete on Mac) key to go back to the main Julia prompt. **Steps 3-5 should only need to be done once per Julia installation.**

## Usage

In order to use the functions, in each new Julia session you must enter

```julia
using FundamentalsNumericalComputation
```

None of the functions are exported, so they must all be prefixed with the package name. However, the constant `FNC` is set as an alias to the package name to make typing more reasonable, e.g., `FNC.lufact`, `FNC.rk23`, etc. 

There is a bare-bones [documentation site][docs-stable-url], but it only gives summaries of the documentation strings that are found in the text. The textbook is meant to be the real guide.

## Startup speed

After installation, importing the package with a `using` statement should only take 5-10 seconds in Julia 1.6 or later. The first plot created in a Julia session may take 20-60 seconds to appear. This is a well-known irritant in the Julia ecosystem, and progress is steadily being made toward reducing the lag. In the meantime, there are a few options available for power users.

### PackageCompiler

The [`PackageCompiler`](https://julialang.github.io/PackageCompiler.jl/) package allows you to compile a new version of the Julia binary to get better startup performance. It's not too hard to use, but it helps if you are comfortable with command line/terminal interfaces.

### Minimal package version

You can speed up loading time by installing a special branch of this package. At the prompt, type the `]` character to enter package mode. If you already installed the default version, enter

```julia
rm FundamentalsNumericalComputation
```

Then enter

```julia
add https://github.com/fncbook/FundamentalsNumericalComputation.jl#fast-load
```

Now, each import of this package should take ten seconds or less. However, you will need to manually install and then load packages in order to run the demos and solve many of the exercises in the book. (This is how Julia is normally used in practice.) In particular, only `LinearAlgebra`, `Polynomials`, `SparseArrays`, and a few packages for loading files and displaying output are provided. Other packages used in the text, and their key dependents, are given in the following table.

| Package name            | Key dependent functions           |
| :---------------------- | :-------------------------------- |
| `Arpack`                | `eigs`                            |
| `Dierckx`               | `Splline1D`                       |
| `DifferentialEquations` | `solve`                           |
| `FFTW`                  | `fft`                             |
| `GraphRecipes`          | `graphplot`                       |
| `Images`                | image loading and display, `Gray` |
| `IncompleteLU`          | `iLU`                             |
| `IterativeSolvers`      | `gmres`, `minres`, `cg`           |
| `LinearMaps`            | `LinearMap`                       |
| `MatrixDepot`           | `matrixdepot`                     |
| `NLsolve`               | `nlsolve`                         |
| `Plots`                 | `plot`, `scatter`, `contour`      |
| `Preconditioners`       | `DiagonalPreconditioner`          |
| `QuadGK`                | `quadgk`                          |
| `SpecialFunctions`      | `besselj`, `gamma`                |
| `TestImages`            | `testimage`                       |

## Alternative numerical packages

These codes are for instructional purposes. They are not recommended for applications. Compared to superior alternatives, they lack generality, efficiency, and robustness to syntax and more subtle mistakes. My personal recommendations for preferable alternatives are as follows. Most of them are demonstrated in the textbook.

| Problem type                                                                                            | Associated Julia packages                                                                                                                                                                                                                                                                         |
| ------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Linear system / least squares                                                                           | [`LinearAlgebra`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#man-linalg) (standard library)                                                                                                                                                                                           |
| Sparse matrix                                                                                           | [SparseArrays](https://docs.julialang.org/en/v1/stdlib/SparseArrays/#Sparse-Arrays), [IterativeSolvers](https://iterativesolvers.julialinearalgebra.org/stable/),  [Arpack](https://arpack.julialinearalgebra.org/stable/), [Preconditioners](https://github.com/mohamed82008/Preconditioners.jl) |
| Polynomial interpolation/approximation                                                                  | [Polynomials](https://juliamath.github.io/Polynomials.jl/stable/), [ApproxFun](https://juliaapproximation.github.io/ApproxFun.jl/stable/)                                                                                                                                                         |
| Polynomial roots                                                                                        | [Polynomials](https://juliamath.github.io/Polynomials.jl/stable/#Root-finding-1)                                                                                                                                                                                                                  |
| Rootfinding                                                                                             | [NLsolve](https://github.com/JuliaNLSolvers/NLsolve.jl)                                                                                                                                                                                                                                           |
| Finite differences                                                                                      | [FiniteDifferences](https://juliadiff.org/FiniteDifferences.jl/latest/), [FiniteDiff](https://github.com/JuliaDiff/FiniteDiff.jl)                                                                                                                                                                 |
| Integration                                                                                             | [QuadGK](https://juliamath.github.io/QuadGK.jl/stable/)                                                                                                                                                                                                                                           |
| Spline                                                                                                  | [Dierckx](https://github.com/kbarbary/Dierckx.jl)                                                                                                                                                                                                                                                 |
| [Initial-value problem](https://diffeq.sciml.ai/latest/tutorials/ode_example/#ode_example)              | [DifferentialEquations](https://diffeq.sciml.ai/latest/)                                                                                                                                                                                                                                          |
| [Boundary-value problem](https://diffeq.sciml.ai/latest/tutorials/bvp_example/#Boundary-Value-Problems) | [DifferentialEquations](https://diffeq.sciml.ai/latest/)                                                                                                                                                                                                                                          |
| Method of lines                                                                                         | [MethodOfLines](https://methodoflines.sciml.ai/dev/) (still under development)                                                                                                                                                                                                                    |


## License

This code stored on this site is under an MIT license. Please see the LICENSE file for details.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://fncbook.github.io/FundamentalsNumericalComputation.jl/stable
