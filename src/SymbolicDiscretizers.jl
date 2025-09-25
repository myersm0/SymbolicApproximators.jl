
module SymbolicDiscretizers

using LinearAlgebra
using Distributions
using StatsBase

include("discretizers.jl")
export SymbolicDiscretizer, alphabet, breakpoints, cardinality

include("utilities.jl")
export numerosity_reduction

include("words.jl")
export Word

# main algorithms
include("sax.jl")
include("ordinal.jl")
include("delta.jl")
export SAX, OrdinalDiscretizer, DeltaDiscretizer
export discretize, reconstruct, distance

include("windows.jl")
export StreamingDiscretizer, update!, get_symbols

end
