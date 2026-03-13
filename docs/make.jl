using Documenter
using PenguinTransportDiffusion

makedocs(
    modules = [PenguinTransportDiffusion],
    authors = "PenguinxCutCell contributors",
    sitename = "PenguinTransportDiffusion.jl",
    format = Documenter.HTML(
        canonical = "https://PenguinxCutCell.github.io/PenguinTransportDiffusion.jl",
        repolink = "https://github.com/PenguinxCutCell/PenguinTransportDiffusion.jl",
        collapselevel = 2,
    ),
    doctest = false,
    pages = [
        "Home" => "index.md",
        "Advection-Diffusion Model" => "transport_diffusion.md",
        "Moving Geometry" => "moving.md",
        "Algorithms" => "algorithms.md",
        "API" => "api.md",
        "Examples" => "examples.md",
    ],
    pagesonly = true,
    warnonly = false,
    remotes = nothing,
)

if get(ENV, "CI", "") == "true"
    deploydocs(
        repo = "github.com/PenguinxCutCell/PenguinTransportDiffusion.jl",
        push_preview = true,
    )
end
