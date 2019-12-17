using Documenter, BioXP

makedocs(
    modules = [BioXP],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Harrison Smith",
    sitename = "BioXP.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/hbsmith/BioXP.jl.git",
)
