push!(LOAD_PATH,"../src/")

using Documenter
using EpiSiming

makedocs(
    sitename = "EpiSiming.jl",
    pages = [
        "Home" => "index.md",
        "types.md",
        "evolution.md",
        "examples.md",
        "api.md"
    ],
    # The following `format` directive is for local browsing, to properly resolve the links
    # Check out the last note in [Building an Empty Document](https://juliadocs.github.io/Documenter.jl/stable/man/guide/#Building-an-Empty-Document)
    format = Documenter.HTML(prettyurls = false),
    modules = [EpiSiming]
)