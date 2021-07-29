module FundamentalsNumericalComputation

export FNC
FNC = FundamentalsNumericalComputation

using Reexport

@info "Re-exporting multiple packages..."

@reexport using LinearAlgebra
@reexport using Statistics
@reexport using SparseArrays
@reexport using Polynomials
@reexport using NLsolve
@reexport using DifferentialEquations
@reexport using Plots
@reexport using PrettyTables
@reexport using Dierckx
@reexport using QuadGK 
@reexport using MatrixDepot
@reexport using GraphRecipes
@reexport using Images
@reexport using TestImages
@reexport using Arpack
@reexport using IterativeSolvers
@reexport using LinearMaps
@reexport using IncompleteLU
@reexport using Preconditioners
@reexport using FFTW
@reexport using FileIO
@reexport using SpecialFunctions
@reexport using LaTeXStrings
@reexport using Printf
@reexport using JLD2
@reexport using Printf

include("chapter01.jl")
include("chapter02.jl")
include("chapter03.jl")
include("chapter04.jl")
include("chapter05.jl")
include("chapter06.jl")
include("chapter08.jl")
include("chapter09.jl")
include("chapter10.jl")
include("chapter11.jl")
include("chapter13.jl")

end
