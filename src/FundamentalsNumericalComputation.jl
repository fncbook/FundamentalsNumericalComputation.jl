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
@reexport using Dierckx: Spline1D
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

"""
This function sets up graphics and other settings used to process the files to 
construct the textbook.
"""
function init_format(;backend=:GR)
    if backend==:GR
        gr()
        default(
            titlefont=(11,"Helvetica"),
            guidefont=(11,"Helvetica"),
            linewidth = 2,
            markersize = 3,
            msa = 0,
            size=(500,320),
            label="",
            html_output_format = "svg"
        )
    elseif backend==:pyplot 
        pyplot()
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.sans-serif"] = append!(["Helvetica"],rcParams["font.sans-serif"])
        rcParams["lines.linewidth"] = 2
        rcParams["lines.markersize"] = 4
        default(
            msa=0,
            size=(500,320),
            label="",
            html_output_format = "svg"
            )
    end
end

end  # module
