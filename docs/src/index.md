# Fundamentals of Numerical Computation

These are the Julia codes from the book [*Fundamentals of Numerical Computation*](https://fncbook.github.io/fnc) by T. A. Driscoll and R. J. Braun.

## Installation

1. [Install Julia.](https://julialang.org/downloads/) Any of the available methods should be fine. The Julia Pro version comes with an integrated editor and many preinstalled packages, and it might be the best choice for those not comfortable with their command line shell.
2. Start Julia on your machine. Look for the `julia>` prompt.
3. At the prompt type the `]` character. This will turn the prompt a different color, indicating that you can give commands to the package manager.
4. At the newly colored prompt, type

   ```julia
   add FundamentalsNumericalComputation
   ```

   The process will take several minutes.
5. Hit the backspace (delete on Mac) key to go back to the main Julia prompt. **Steps 3-5 should only need to be done once per installation.**

## Usage

In order to use the functions, in each new Julia session you must enter

```julia
using FundamentalsNumericalComputation
```

After this completes, all the text's functions can be accessed with the prefix `FNC`. E.g., `FNC.lufact`, `FNC.rk23`, etc.

