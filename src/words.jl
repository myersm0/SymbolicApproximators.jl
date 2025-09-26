
struct Word{A, T} <: AbstractVector{T}
	approximator::Ref{A}
	symbols::Vector{T}
end

function Word(ca::CA, values) where CA <: ContinuousApproximator
	T = eltype(ca)
	return Word{CA, T}(Ref(ca), values)
end

function Word(sa::SA, values) where SA <: SymbolicApproximator
	α = sa.α
	β = sa.β
	T = eltype(α)
	n = length(values)
	symbols = Vector{T}(undef, n)
	@inbounds @simd for i in eachindex(values)
		v = Float64(values[i])
		index = searchsortedlast(β, v)
		symbols[i] = α[index + 1]
	end
	return Word{SA, T}(Ref(sa), symbols)
end

alphabet(word::Word) = alphabet(word.approximator[])
breakpoints(word::Word) = breakpoints(word.approximator[])
cardinality(word::Word) = cardinality(word.approximator[])

Base.size(w::Word) = size(w.symbols)
Base.getindex(word::Word, args...) = getindex(word.symbols, args...)
Base.IndexStyle(::Type{<:Word}) = IndexLinear()
Base.iterate(w::Word, state = 1) = iterate(w.symbols, state)





