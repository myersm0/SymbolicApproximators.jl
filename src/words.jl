
struct Word{A <: AbstractApproximator, T}
	approximator::Ref{A}
	data::Vector{T}
end

function Word(ca::CA, values::AbstractVector) where CA <: ContinuousApproximator
	return Word{CA, Float64}(Ref(ca), float_values)
end

function Word(sa::SA, values::AbstractVector) where SA <: SymbolicApproximator
	β = breakpoints(sa)
	n = length(values)
	indices = Vector{Int}(undef, n)
	@inbounds @simd for i in eachindex(values)
		v = Float64(values[i])
		indices[i] = searchsortedlast(β, v)
	end
	return Word{SA, Int}(Ref(sa), indices)
end

Base.keys(w::Word{<:ContinuousApproximator, Float64}) = Base.OneTo(length(w.data))
Base.keys(w::Word{<:SymbolicApproximator, Int}) = w.data

function Base.values(w::Word{<:SymbolicApproximator, Int})
	α = alphabet(w)
	T = eltype(α)
	n = length(w.data)
	symbols = Vector{T}(undef, n)
	@inbounds @simd for i in 1:n
		symbols[i] = α[w.data[i]]
	end
	return symbols
end

Base.values(w::Word{<:ContinuousApproximator, Float64}) = w.data

alphabet(w::Word) = alphabet(w.approximator[])
breakpoints(w::Word) = breakpoints(w.approximator[])
cardinality(w::Word) = cardinality(w.approximator[])
word_size(w::Word) = word_size(w.approximator[])

Base.length(w::Word) = length(w.data)
Base.size(w::Word) = (length(w),)
Base.eltype(::Word{A, T}) where {A, T} = T

function Base.getindex(w::Word{<:ContinuousApproximator}, args...)
	return getindex(w.data, args...)
end

function Base.getindex(w::Word{<:SymbolicApproximator}, args...)
	@boundscheck checkbounds(w.data, args...)
	return getindex(alphabet(w), getindex(w.data, args...))
end

function Base.iterate(w::Word{<:ContinuousApproximator}, state = 1)
	state > length(w) && return nothing
	return (w.data[state], state + 1)
end

function Base.iterate(w::Word{<:SymbolicApproximator}, state = 1)
	state > length(w) && return nothing
	return (alphabet(w)[w.data[state]], state + 1)
end

Base.lastindex(w::Word) = length(w)

Base.:(==)(w1::Word, w2::Word) = 
	w1.approximator[] == w2.approximator[] && w1.data == w2.data

Base.hash(w::Word, h::UInt) = hash(w.data, hash(w.approximator[], h))

