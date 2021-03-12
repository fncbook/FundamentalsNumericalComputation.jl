# Fundamentals of Numerical Computation

These are source files used to generate the online text [*Fundamentals of Numerical Computation*](https://fncbook.github.io/fnc) by T. A. Driscoll and R. J. Braun.

## Using the book in a course

It is advisable to use one of rows of this table for stable URLs whose content will not change except for typos and essential fixes.

|  Version  |  Date created  |  Text website  |  Source files  |
|-----------|----------------|----------------|----------------|
| v1.0  | July 2020 | [fncbook.github.io/v1.0](https://fncbook.github.io/v1.0) | [zip and tar.gz](https://github.com/fncbook/fnc/releases/tag/v1.0) |

## Installing files locally

Each chapter includes notebooks that were run with Julia 1.4. They should work with any later 1.x release as well (provided dependencies honor the standard versioning scheme).

To get the book functions as well as all demos:

1. [Install Julia.](https://julialang.org/downloads/) Any of the available methods should be fine. The Julia Pro version comes with an integrated editor and many preinstalled packages, and it might be the best choice for those not comfortable with their command line shell.
2. Start Julia on your machine. Look for the `julia>` prompt.
3. At the prompt type the `]` character. This will turn the prompt a different color, indicating that you can give commands to the package manager.
4. At the newly colored prompt, type

   ```julia
   add https://github.com/fncbook/fnc
   ```

   The process will take several minutes.
5. Hit the backspace (delete on Mac) key to go back to the main Julia prompt. **Steps 3-5 should only need to be done once per installation.**
6. In order to use the functions, in each new Julia session you must enter

   ```julia
   using FundamentalsNumericalComputation
   ```

   After this completes, all the text's functions can be accessed with the prefix `FNC`. E.g., `FNC.lufact`, `FNC.rk23`, etc.

### Startup speed

There are two sources of noticeable slowness:

1. Any `using` statement requires a compilation step after a package or its dependencies have been updated. As of July 2020, this can consume several or even tens of minutes for `FNC`, due mainly to its dependence on [`Plots`](http://docs.juliaplots.org/latest/) and [`DifferentialEquations`](https://docs.sciml.ai/latest/index.html). Even on subsequent invocations, `using` this package can take 10-30 seconds to load.
2. The first plot created in a Julia session may take up to a minute, which is typical for `Plots`.

Both of these lag issues are well known to Julia developers and are the target of future improvements. In the meantime power users might investigate using [`PackageCompiler`](https://julialang.github.io/PackageCompiler.jl/dev/) to get vastly better performance on both issues.

## Tooling

The book was created using [jupyter-book](https://jupyterbook.org) and [VS Code](https://code.visualstudio.com/) with the [Julia extension](https://github.com/julia-vscode/julia-vscode). These tools are already excellent and continue to evolve rapidly.

## License

All the material is copyrighted. The text is derived from *Fundamentals of Numerical Computation* by T. A. Driscoll and R. J. Braun, copyright 2017 by the Society for Industrial and Applied Mathematics, with all rights reserved. The Julia exercises are licensed under CC Attribution-ShareAlike 4.0. The code generating the site is under an MIT license. Please see the LICENSE file for details.
