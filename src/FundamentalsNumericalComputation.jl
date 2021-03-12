module FundamentalsNumericalComputation

export FNC
FNC = FundamentalsNumericalComputation

using Reexport

@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using Polynomials
@reexport using NLsolve
@reexport using OrdinaryDiffEq
@reexport using Plots
@reexport using PrettyTables
@reexport using QuadGK 
@reexport using Interpolations
@reexport using GraphRecipes
@reexport using Images
@reexport using TestImages
@reexport using Arpack
@reexport using IterativeSolvers
@reexport using LinearMaps
@reexport using IncompleteLU
@reexport using Preconditioners
@reexport using OffsetArrays
@reexport using FFTW
@reexport using FileIO
@reexport using SpecialFunctions
@reexport using LaTeXStrings
@reexport using JLD2


@info "Re-exporting multiple packages"

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

function magic(n)
	# Adapted from an implementation by Phillip Berndt
	if isodd(n)
		p = 1:n
		A = p .- div(n+3, 2) .+ p'
		B = 2p .- 2 .+ p'
		return n*mod.(A,n) + mod.(B,n) .+ 1
	elseif iseven(n) && iseven(div(n,2))
		J = @. div((1:n)%4, 2)
		K = J' .== J
		M = (0:n-1) .+ (1:n:n^2)'
		M[K] = n^2 + 1 .- M[K]
		return M
	else
		p = div(n,2)
		M = magic(p)
		M = [M M.+2p^2; M.+3p^2 M.+p^2]
		(n == 2) && return M
		i = 1:p
		k = div(n-2,4)
		j = [(1:k); ((n-k+2):n)]
		M[[i; i.+p],j] = M[[i.+p; i],j]
		i = k+1
		j = [1; i]
		M[[i; i+p],j] = M[[i+p; i],j]
		return M
	end
end

end
