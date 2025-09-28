
module SymbolicApproximators

using Distances
using Distributions
using StatsBase
using StaticArrays
using OffsetArrays

include("approximators.jl")
export Approximator, ContinuousApproximator, SymbolicApproximator
export alphabet, breakpoints, cardinality, alphabet_size, word_size

include("words.jl")
export Word, width, compression_rate

include("encode.jl")
export encode

# main algorithms
include("paa.jl")
include("sax.jl")
include("esax.jl")
#include("ordinal.jl")
#include("delta.jl")
export PAA, SAX, ESAX

include("distances.jl")
export evaluate

include("show.jl")

end

