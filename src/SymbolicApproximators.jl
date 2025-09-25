
module SymbolicApproximators

using LinearAlgebra
using Distributions
using StatsBase

include("approximators.jl")
export Approximator, ContinuousApproximator, SymbolicApproximator
export alphabet, breakpoints, cardinality, alphabet_size

include("words.jl")
export Word

# main algorithms
include("paa.jl")
include("sax.jl")
include("ordinal.jl")
include("delta.jl")
export PAA, SAX, OrdinalApproximator, DeltaApproximator
export discretize, reconstruct, distance

include("windows.jl")
export StreamingApproximator, update!, get_symbols

end
