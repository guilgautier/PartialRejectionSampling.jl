# Sampling methods

Below is the list of point processes which can be sampled exactly using `PartialRejectionSampling.jl`

* [Spatial point processes](@ref)
  * [`PRS.HardCorePointProcess`](@ref)
  * [`PRS.HomogeneousPoissonPointProcess`](@ref)
  * [`PRS.StraussPointProcess`](@ref)
* [Graph point processes](@ref)
  * [`PRS.HardCoreGraph`](@ref)
  * [`PRS.Ising`](@ref)
  * [`PRS.RootedSpanningForest`](@ref)
  * [`PRS.SinkFreeGraph`](@ref)
* Miscellaneous
  * [`PRS.PatternFreeString`](@ref)

## Default exact sampler

```@docs
PRS.generate_sample
```

## Partial Rejection Sampling (PRS)

[GuJeLi19](@cite) developed the PRS methodology to generate exact samples from product distributions of the form ``\otimes_{n=1}^{N} \mu_n`` subject to some constraints.
It requires access to an exact sampler for each ``\mu_n`` and an oracle to check the violation of the constraints.

Given an initial sample from ``\otimes_{n=1}^{N} \mu_n``:

* [Vanilla rejection sampling](https://en.wikipedia.org/wiki/Rejection_sampling) resample all variables if any constraint is violated; until all constraints are satisfied,
* Partial rejection sampling instead constructs a subset of variables to be resampled, starting from variables involved in violated constraints, and preserves the state of the variables outside of this resampling set; until all constraints are satisfied.

In both cases, the output sample is guaranteed to have the right distribution, i.e., the product distribution subject to the prescribed constraints.

```@docs
PRS.generate_sample_prs
```

## Grid Partial Rejection Sampling (grid PRS)

[MoKr20](@cite) adapted the idea of [Partial Rejection Sampling (PRS)](@ref) (originally derived in the finite setting) to generate exact samples from [Spatial point processes](@ref) with finite range of interaction, noted ``r > 0``.

The name of the method *grid* PRS results from the combination of the [Partial Rejection Sampling (PRS)](@ref) methodology with a partitioning of the domain of ``\mathbb{R}^2`` where the target [`PRS.AbstractSpatialPointProcess`](@ref) -- is divided into ``N`` cells of type [`PRS.SpatialCellGridPRS`](@ref) with size ``r \times r``,

To draw the correspondence with framework of PRS developed by [GuJeLi19](@cite), ``\mu_n`` represents the target [`PRS.AbstractSpatialPointProcess`](@ref) restricted to the ``n``-th cell, i.e., variables are point configurations.

**Note** that the grid PRS methodology requires an efficient sampling algorithm to generate exact samples on each cell involved in the partitioning of original domain, see, e.g., [`PRS.generate_sample_dcftp`](@ref) and [`PRS.generate_sample_prs`](@ref).

**See also** closely related variants of grid PRS for

* Spatial point processes
  * [Hub20](@cite)
* Graphical models
  * [Ising model](@ref)
  * [FeViYi19](@cite), [FeGuYi19](@cite)

```@autodocs
Modules = [PartialRejectionSampling]
Pages   = ["grid_prs.jl"]
Private = true
Order = [:module, :constant, :type, :function, :macro]
```

## Uniform sampling in spatial windows

```@docs
Base.rand
```

## Dominated Coupling From The Past (dCFTP)

Implementation of dominated Coupling From The Past (dCFTP) developed by [KeMo99](@cite) and [KeMo00](@cite)

**See also**

* [Hub16](@cite)
* [Kendall's notes on perfect simulation](https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/kendall/personal/ppt/428.pdf)

```@autodocs
Modules = [PartialRejectionSampling]
Pages   = ["dominated_cftp.jl"]
Private = true
Order = [:module, :constant, :type, :function, :macro]
```
