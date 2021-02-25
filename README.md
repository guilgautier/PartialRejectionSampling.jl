# PartialRejectionSampling.jl

<!-- [![][docs-stable-img]][docs-stable-url]  -->
[![][docs-dev-img]][docs-dev-url]

This module provides a [Julia](https://julialang.org/) implementation of the Partial Rejection Sampling (PRS) methodology, recently developed by [Guo, Jerrum and Liu (2019)](https://guilgautier.github.io/PartialRejectionSampling.jl/dev/references/).
With PRS, you generate exact samples from product distributions subject to some constraints, see e.g., some [Graph point processes](@ref) and [Spatial point processes](@ref).

Given an initial sample from the (unconstrained) product distribution:

- [Vanilla rejection sampling](https://en.wikipedia.org/wiki/Rejection_sampling) resample all variables if any constraint is violated; until all constraints are satisfied,
- Partial rejection sampling instead identifies a subset of variables to be resampled, starting from variables involved in violated constraints, and preserves the state of the variables outside of this resampling set; until all constraints are satisfied.

In both cases, the output sample is guaranteed to have the right distribution, i.e., the product distribution subject to the prescribed constraints.

## Getting Started

### Installation

[`PartialRejectionSampling.jl`](https://github.com/guilgautier/PartialRejectionSampling.jl) is not a registered package, yet.
Nevertheless, you can to install it through

```julia
julia> ]add https://github.com/guilgautier/PartialRejectionSampling.jl
```

see also [how to manage packages with `Pkg`](https://julialang.github.io/Pkg.jl/stable/managing-packages/##Adding-packages-1).

### Usage

To start using the package, simply enter

```julia
julia> using PartialRejectionSampling
# const PRS = PartialRejectionSampling is made available so you can then use
# PRS.<type/function_you_want_to_use>
```

### Tutorial Jupyter notebooks

You can have a look at the tutorial Jupyter [`notebooks`](https://github.com/guilgautier/PartialRejectionSampling.jl/blob/master/notebooks) to play with code.

### Documentation

The documentation is currently available at [![][docs-dev-img]][docs-dev-url]

## Bug reports - Contributions

Feel free to [raise issues](https://github.com/guilgautier/PartialRejectionSampling.jl/issues), make comments or make [pull requests](https://github.com/guilgautier/PartialRejectionSampling.jl/pulls).
Any feedback is welcome :smiley:.

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://guilgautier.github.io/PartialRejectionSampling.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://guilgautier.github.io/PartialRejectionSampling.jl/stable
