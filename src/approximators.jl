
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

function Base.isequal(sa1::SymbolicApproximator, sa2::SymbolicApproximator)
	return false
end

function Base.isequal(sa1::T, sa2::T) where T <: SymbolicApproximator
	return word_size(sa1) == word_size(sa2) && alphabet_size(sa1) == alphabet_size(sa2)
end

function Base.isequal(ca1::T, ca2::T) where T <: ContinuousApproximator
	return word_size(ca1) == word_size(ca2)
end

Base.:(==)(sa1::SymbolicApproximator, sa2::SymbolicApproximator) =
	isequal(sa1, sa2)



