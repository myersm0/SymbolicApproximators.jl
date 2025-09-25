
module SymbolicApproximators

using LinearAlgebra
using Distributions
using StatsBase

include("approximators.jl")
export SymbolicApproximator, alphabet, breakpoints, cardinality

include("utilities.jl")
export numerosity_reduction

include("words.jl")
export Word

# main algorithms
include("sax.jl")
include("ordinal.jl")
include("delta.jl")
export SAX, OrdinalApproximator, DeltaApproximator
export discretize, reconstruct, distance

include("windows.jl")
export StreamingApproximator, update!, get_symbols

end
