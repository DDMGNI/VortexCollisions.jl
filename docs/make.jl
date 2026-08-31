using VortexCollisions
using Documenter

DocMeta.setdocmeta!(VortexCollisions, :DocTestSetup, :(using VortexCollisions); recursive = true)

makedocs(;
    modules = [VortexCollisions],
    authors = "Michael Kraus",
    repo = "https://github.com/DDMGNI/VortexCollisions.jl/blob/{commit}{path}#{line}",
    sitename = "VortexCollisions.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://DDMGNI.github.io/VortexCollisions.jl",
        edit_link = "master",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Library" => "library.md"
    ]
)

deploydocs(;
    repo = "github.com/DDMGNI/VortexCollisions.jl",
    devbranch = "master"
)
