using Documenter
using FundamentalsNumericalComputation

makedocs(
    sitename = "FundamentalsNumericalComputation",
    format = Documenter.HTML(),
    modules = [FundamentalsNumericalComputation],
    pages = [
        "index.md",
        "Functions" => "functions.md",
 #       "Square linear systems" => "chapter02.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
