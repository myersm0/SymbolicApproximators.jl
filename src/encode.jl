
function windows(
		values::AbstractVector, window_size::Int; 
		stride = 1, delay = 1
	)
	n = length(values)
	max_start = n - delay * (window_size - 1)
	return (
		@view(
			values[i:delay:i+delay*(window_size-1)]
		)
		for i in 1:stride:max_start
	)
end

function segments(values::AbstractVector, n_segments::Int)
	n = length(values)
	segment_length = n ÷ n_segments
	return (
		@view(
			values[(i-1)*segment_length+1 : (i==n_segments ? n : i*segment_length)]
		)
		for i in 1:n_segments
	)
end

function segments(values::AbstractVector, segment_lengths::Vector{Int})
	n = length(values)
	sum(segment_lengths) == n || 
		error("Segment lengths must sum to values length $n")
	start = 1
	return map(segment_lengths) do len
		segment = @view(values[start:start+len-1])
		start += len
		segment
	end
end

"""
    encode(approximator, values)

Turn continuous `values` into a lower-dimensional approximation.

The configuration defined in `approximator::AbstractApproximator` will
determine the algorithm and parameters for doing so.

Returns a `Word`.
"""
function encode(a::AbstractApproximator, values::AbstractVector)
	n = length(values)
	w = word_size(a)
	segs = segments(values, w)
	return Word(a, [_encode_segment(a, seg) for seg in segs], n)
end

function encode!(
		a::AbstractApproximator, dest::C, values::AbstractVector
	) where C <: AbstractArray
	n = length(values)
	w = word_size(sa)
	segs = segments(values, w)
	return Word(
		a, dest, [_encode_segment(a, seg) for seg in segs], n
	)
end

approximate(args...) = encode(args...)

approximate!(args...) = encode!(args...)




