getRNG(seed::Integer = -1) = seed >= 0 ? Random.MersenneTwister(seed) : Random.GLOBAL_RNG
getRNG(seed::Union{Random.MersenneTwister,Random._GLOBAL_RNG}) = seed

sigmoid(x) = @. 1 / (1 + exp(-x))
