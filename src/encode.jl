
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

function encode(sa::ContinuousApproximator, values::AbstractVector)
	w = word_size(sa)
	segs = segments(values, w)
	return [_encode_segment(sa, seg) for seg in segs]
end

function encode(sa::SymbolicApproximator, values::AbstractVector)
	w = word_size(sa)
	segs = segments(values, w)
	return Word(sa, [_encode_segment(sa, seg) for seg in segs])
end






