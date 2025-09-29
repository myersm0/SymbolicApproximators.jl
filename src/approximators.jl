
abstract type AbstractApproximator end
abstract type ContinuousApproximator <: AbstractApproximator end
abstract type SymbolicApproximator <: AbstractApproximator end

word_size(a::AbstractApproximator) = a.w
alphabet(sa::ContinuousApproximator) = Inf
breakpoints(sa::ContinuousApproximator) = nothing
cardinality(sa::ContinuousApproximator) = nothing
alphabet_size(sa::ContinuousApproximator) = Inf

alphabet(sa::SymbolicApproximator) = range('a'; length = sa.a)
breakpoints(sa::SymbolicApproximator) = sa.β
cardinality(sa::SymbolicApproximator) = sa.a
alphabet_size(sa::SymbolicApproximator) = cardinality(sa)

"""
    encode(approximator, values)

Approximate continuous `values` in lower-dimensional representation.

The configuration defined in `approximator::AbstractApproximator` will
determine the algorithm and parameters for doing so.

Returns a `Word`.
"""
function encode(::AbstractApproximator, values::AbstractVector) end

approximate(sa::AbstractApproximator, values) = encode(sa, values)




