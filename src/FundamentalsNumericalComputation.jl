module FundamentalsNumericalComputation

export FNC
FNC = FundamentalsNumericalComputation

using Reexport

# Required by the package itself
@reexport using Polynomials
@reexport using OrdinaryDiffEq
@reexport using LinearAlgebra
@reexport using SparseArrays

# Not strictly required, but among the more obscure imports in the text, and fast to load.
@reexport using FileIO
@reexport using LaTeXStrings
@reexport using JLD2
@reexport using Printf
@reexport using PrettyTables

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
