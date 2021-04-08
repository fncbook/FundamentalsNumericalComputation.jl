using Documenter
using FundamentalsNumericalComputation

Documenter.Writers.HTMLWriter.HTML(sidebar_sitename=false)
makedocs(
    sitename = "FNC Functions",
    format = Documenter.HTML(),
    modules = [FundamentalsNumericalComputation],
    pages = [
        "index.md",
        "Functions" => "functions.md",
    ]
)

deploydocs(
    repo = "github.com/fncbook/FundamentalsNumericalComputation.jl.git",
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
