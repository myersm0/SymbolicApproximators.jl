
module SymbolicApproximators

using Distances
using Distributions
using StatsBase

include("approximators.jl")
export Approximator, ContinuousApproximator, SymbolicApproximator
export alphabet, breakpoints, cardinality, alphabet_size

include("words.jl")
export Word

include("encode.jl")
export encode

# main algorithms
include("paa.jl")
include("sax.jl")
#include("ordinal.jl")
#include("delta.jl")
export PAA, SAX #, OrdinalApproximator, DeltaApproximator

include("distances.jl")
export evaluate

include("show.jl")

end

