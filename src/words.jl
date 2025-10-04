
abstract type WordStyle end
struct SimpleWord <: WordStyle end
struct MultiWord{Int} <: WordStyle end

width(::SimpleWord) = 1
width(::MultiWord{W}) where W = W

struct Word{A <: AbstractApproximator, T <: AbstractArray, W}
	approximator::Ref{A}       # ref to what generated the Word
	data::T                    # the main content (varies by approximator)
	n::Int                     # number of timepoints in the original time series
end

function Word(a::AbstractApproximator, args...)
	return Word(WordStyle(a), a, args...)
end

function Word(
		::SimpleWord, ca::CA, values::AbstractVector, n::Integer
	) where CA <: ContinuousApproximator
	return Word{CA, Vector{Float64}, 1}(Ref(ca), values, n)
end

function Word(
		::SimpleWord, sa::SA, values::AbstractVector, n::Integer
	) where SA <: SymbolicApproximator
	β = breakpoints(sa)
	indices = Vector{Int}(undef, length(values))
	@inbounds @simd for i in eachindex(values)
		v = Float64(values[i])
		indices[i] = searchsortedlast(β, v)
	end
	return Word{SA, Vector{Int}, 1}(Ref(sa), indices, n)
end

function Word(
		::SimpleWord, sa::SA, dest::AbstractArray, values::AbstractVector, n::Integer
	) where SA <: SymbolicApproximator
	β = breakpoints(sa)
	@inbounds @simd for i in eachindex(values)
		dest[i] = searchsortedlast(β, Float64(values[i]))
	end
	return Word{SA, typeof(dest), 1}(Ref(sa), dest, n)
end

function Word(
		::MultiWord{W}, sa::SA, values, n::Integer
	) where {W, SA <: SymbolicApproximator}
	β = breakpoints(sa)
	indices = Vector{SVector{W, Int}}(undef, length(values))
	for (i, v) in enumerate(values)
		indices[i] = SVector{W, Int}(searchsortedlast(β, vi) for vi in v)
	end
	return Word{SA, Vector{SVector{W, Int}}, W}(Ref(sa), indices, n)
end

WordStyle(w::Word{A, T, 1}) where {A, T} = SimpleWord()
WordStyle(w::Word{A, T, W}) where {A, T, W} = MultiWord{W}()

Base.keys(w::Word{<:ContinuousApproximator}) = Base.OneTo(length(w.data))
Base.keys(w::Word{<:SymbolicApproximator}) = w.data

"""
    values(w::Word)

Get a vector of w's symbolic representation.
"""
function Base.values(w::Word)
	return values(WordStyle(w), w)
end
	
function Base.values(::SimpleWord, w::Word{<:SymbolicApproximator})
	return 'a' .+ w.data
end

function Base.values(::MultiWord{W}, w::Word{<:SymbolicApproximator}) where W
	return ['a' .+ v for v in w.data]
end

Base.values(w::Word{<:ContinuousApproximator}) = w.data

alphabet(w::Word) = alphabet(w.approximator[])
breakpoints(w::Word) = breakpoints(w.approximator[])
cardinality(w::Word) = cardinality(w.approximator[])
alphabet_size(w::Word) = alphabet_size(w.approximator[])
word_size(w::Word) = word_size(w.approximator[])
width(w::Word{A, T, W}) where {A, T, W} = W
compression_rate(n, w, width = 1) = n / (w * width)
compression_rate(w::Word) = compression_rate(w.n, word_size(w), width(w))

Base.length(w::Word) = word_size(w)
Base.size(w::Word) = (length(w),)
Base.eltype(::Word{A, T, W}) where {A, T, W} = eltype(T)

function Base.getindex(w::Word{<:ContinuousApproximator}, args...)
	return getindex(w.data, args...)
end

function Base.getindex(w::Word{<:SymbolicApproximator}, args...)
	return 'a' + getindex(w.data, args...)
end

function Base.iterate(w::Word{<:ContinuousApproximator}, state = 1)
	state > length(w) && return nothing
	return (w.data[state], state + 1)
end

function Base.iterate(w::Word{<:SymbolicApproximator}, state = 1)
	state > length(w) && return nothing
	return ('a' + w.data[state], state + 1)
end

Base.lastindex(w::Word) = length(w)

Base.:(==)(w1::Word, w2::Word) = 
	w1.approximator[] == w2.approximator[] && w1.data == w2.data && w1.n == w2.n

Base.hash(w::Word, h::UInt) = hash(w.data, hash(w.approximator[], h))

