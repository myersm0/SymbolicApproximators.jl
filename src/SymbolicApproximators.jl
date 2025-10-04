
module SymbolicApproximators

using Distances
using Distributions
using StatsBase
using StaticArrays

include("approximators.jl")
export Approximator, ContinuousApproximator, SymbolicApproximator
export alphabet, breakpoints, cardinality, alphabet_size, word_size

include("words.jl")
export Word, width, compression_rate

include("encode.jl")
export encode, encode!

include("methods/paa.jl")
include("methods/sax.jl")
include("methods/esax.jl")
export PAA, SAX, ESAX

include("distances.jl")
export evaluate, MinDist, MINDIST, mindist

include("show.jl")

end

