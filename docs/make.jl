using Documenter, PartialRejectionSampling

# DocMeta.setdocmeta!(PartialRejectionSampling, :DocTestSetup, :(using PartialRejectionSampling); recursive=true)

makedocs(
    sitename="My Documentation",
    # doctest=:fix  # = false to skip doctests
)
