# Spatial point processes

Partial rejection sampling (PRS) can be applied to generate exact samples from pairwise Gibbs point processes with finite interaction range.

For more details on spatial point processes please refer to [MoWa04](@cite) and references therein.

```@docs
PRS.AbstractSpatialPointProcess
PRS.window(::PRS.AbstractSpatialPointProcess)
PRS.dimension(::PRS.AbstractSpatialPointProcess)
```

## Poisson point process

```@autodocs
Modules = [PartialRejectionSampling]
Pages   = ["spatial/poisson.jl"]
Private = true
Order = [:module, :constant, :type, :function, :macro]
```

## Hard core point process

```@autodocs
Modules = [PartialRejectionSampling]
Pages   = ["spatial/hard_core.jl"]
Private = true
Order = [:module, :constant, :type, :function, :macro]
```

## Strauss point process

```@autodocs
Modules = [PartialRejectionSampling]
Pages   = ["spatial/strauss.jl"]
Private = true
Order = [:module, :constant, :type, :function, :macro]
```

## Windows

```@docs
PRS.AbstractSpatialWindow
```

```@autodocs
Modules = [PartialRejectionSampling]
Pages   = ["spatial/window.jl"]
Private = true
Order = [:type, :function, :macro]
```
