using Documenter, PartialRejectionSampling
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "prs.bib"))

# DocMeta.setdocmeta!(PartialRejectionSampling, :DocTestSetup, :(using PartialRejectionSampling); recursive=true)

makedocs(
    bib,
    sitename="PartialRejectionSampling.jl",
    modules = [PartialRejectionSampling],
    pages=[
        "Home"                      => "index.md",
        "Sampling methods"          => "sampling.md",
        "Graph point processes"     => "graph.md",
        "Spatial point processes"   => "spatial.md",
        "Misc"                      => "pattern_free_string.md",
        "References"                => "references.md"
        ]
    # doctest=:fix  # = false to skip doctests
)

deploydocs(
    repo="github.com/guilgautier/PartialRejectionSampling.jl.git",
)
