# `]activate .``
# `using Revise, Documenter`
# `;cd docs`
# `include("make.jl")`
#
push!(LOAD_PATH,Base.current_project(@__DIR__))

using Documenter
using EpiSiming

makedocs(
    sitename = "EpiSiming.jl",
    pages = [
        "Home" => "index.md",
        "Structures" => "structures.md",
        "Evolution" => "evolution.md",
        "Examples" => "examples.md",
        "API" => "api.md"
    ],
    # The following `format` directive is for local browsing, to properly resolve the links
    # Check out the last note in [Building an Empty Document](https://juliadocs.github.io/Documenter.jl/stable/man/guide/#Building-an-Empty-Document)
    format = Documenter.HTML(prettyurls = false),
    modules = [EpiSiming]
)

deploydocs(
    repo = "github.com/ModSiming/EpiSiming.jl.git",
)
