
struct PAA{T <: Real} <: ContinuousApproximator{T}
	w::Int  # word size
end

PAA(w::Integer) = PAA{Float64}(w)

function _encode_segment(a::PAA, values::AbstractVector{<:Real})
	return mean(values)
end

function encode(a::PAA, values::AbstractVector{<:Real})
	n = length(values)
	w = word_size(a)
	segment_length = n ÷ w  # todo: handle case of not divisible
	result = Vector{Float64}(undef, w)
	@inbounds for i in 1:w
		start = round(Int, (i - 1) * segment_length + 1)
		stop = min(n, round(Int, i * segment_length))
		segment = view(values, start:stop)
		result[i] = _encode_segment(a, segment)
	end
	return Word(a, result)
end




