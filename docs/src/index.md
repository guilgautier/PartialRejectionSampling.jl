# [`PartialRejectionSampling.jl`](https://github.com/guilgautier/PartialRejectionSampling.jl) documentation

This module provides a [Julia](https://julialang.org/) implementation of the Partial Rejection Sampling (PRS) methodology developed by [GuJeLi19](@cite).
With PRS, you generate exact samples from product distributions subject to some constraints, see e.g., some [Graph point processes](@ref) and [Spatial point processes](@ref).

Given an initial sample from the (unconstrained) product distribution:

- [Vanilla rejection sampling](https://en.wikipedia.org/wiki/Rejection_sampling) resample all variables if any constraint is violated; until all constraints are satisfied,
- Partial rejection sampling instead identifies a subset of variables to be resampled, starting from variables involved in violated constraints, and preserves the state of the variables outside of this resampling set; until all constraints are satisfied.

In both cases, the output sample is guaranteed to have the right distribution, i.e., the product distribution subject to the prescribed constraints.

## Getting Started

### Installation

[`PartialRejectionSampling.jl`](https://github.com/guilgautier/PartialRejectionSampling.jl) is not yet a registered package. But you can to install it through

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

You can also have a look at the tutorial Jupyter [`notebooks`](https://github.com/guilgautier/PartialRejectionSampling.jl/blob/master/notebooks) to play with the code.

## Index

### Types

```@index
Order   = [:type]
```

### Functions

```@index
Order   = [:function]
```
