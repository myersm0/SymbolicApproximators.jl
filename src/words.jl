
struct Word{A <: AbstractApproximator, T, W}
	approximator::Ref{A}       # ref to what generated the Word
	data::Vector{T}            # the main content (varies by approximator)
	n::Int                     # number of timepoints in the original time series
end

function Word(
		ca::CA, values::AbstractVector, n::Integer
	) where CA <: ContinuousApproximator
	return Word{CA, Float64, width(CA)}(Ref(ca), values, 1, n)
end

function Word(
		sa::SA, values::AbstractVector, n::Integer
	) where SA <: SymbolicApproximator
	β = breakpoints(sa)
	indices = Vector{Int}(undef, length(values))
	@inbounds @simd for i in eachindex(values)
		v = Float64(values[i])
		indices[i] = searchsortedlast(β, v)
	end
	return Word{SA, Int, width(sa)}(Ref(sa), indices, n)
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

function Base.values(w::Word{SA, SVector{3, Int}, 3}) where SA <: SymbolicApproximator
	α = alphabet(w)
	T = eltype(α)
	n = length(w.data)
	symbols = Vector{SVector{3, T}}(undef, n)
	for (i, v) in enumerate(w.data)
		symbols[i] = SVector{3, T}(α[vi] for vi in v)
	end
	return symbols
end

Base.values(w::Word{<:ContinuousApproximator, Float64}) = w.data

alphabet(w::Word) = alphabet(w.approximator[])
breakpoints(w::Word) = breakpoints(w.approximator[])
cardinality(w::Word) = cardinality(w.approximator[])
alphabet_size(w::Word) = alphabet_size(w.approximator[])
word_size(w::Word) = word_size(w.approximator[])
width(w::Word{A, T, W}) where {A, T, W} = W
compression_rate(n, w, width = 1) = sqrt(n / (w * width))
compression_rate(w::Word) = compression_rate(w.n, word_size(w), width(w))

Base.length(w::Word) = word_size(w)
Base.size(w::Word) = (length(w),)
Base.eltype(::Word{A, T, W}) where {A, T, W} = T

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
	w1.approximator[] == w2.approximator[] && w1.data == w2.data && w1.n == w2.n

Base.hash(w::Word, h::UInt) = hash(w.data, hash(w.approximator[], h))

