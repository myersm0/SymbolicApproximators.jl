
abstract type Approximator{T} end
abstract type ContinuousApproximator{T} <: Approximator{T} end
abstract type SymbolicApproximator{T, A} <: Approximator{T} end

abstract type PreprocessingStyle end
struct NoPreprocessing <: PreprocessingStyle end
struct Normalize <: PreprocessingStyle end

word_size(a::Approximator) = a.w
Base.eltype(::Approximator{T}) where T = T

alphabet(sa::ContinuousApproximator{T}) where T = T(Inf)
breakpoints(sa::ContinuousApproximator) = nothing
cardinality(sa::ContinuousApproximator) = nothing
alphabet_size(sa::ContinuousApproximator{T}) where T = T(Inf)

alphabet(sa::SymbolicApproximator) = sa.α
breakpoints(sa::SymbolicApproximator) = sa.β
cardinality(sa::SymbolicApproximator) = length(alphabet(sa))
alphabet_size(sa::SymbolicApproximator) = cardinality(sa)

function encode(::Approximator, values::AbstractVector) end

approximate(sa::Approximator, values) = encode(sa, values)




