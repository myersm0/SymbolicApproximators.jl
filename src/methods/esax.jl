
struct ESAX <: SymbolicApproximator
	w::Union{Nothing, Int}  # word size
	a::UInt32               # alphabet size
	β::Vector{Float64}      # breakpoints
end

function ESAX(w::Integer, a::Integer)
	β = quantile.(Normal(), (1:a) ./ a)
	return ESAX(w, a, β)
end

function _encode_segment(sa::ESAX, values::AbstractVector{<:Real})
	# single pass to find min, max, and sum
	min_value = values[1]
	max_value = values[1]
	min_index = 1
	max_index = 1
	sum_values = values[1]
	@inbounds for i in 2:length(values)
		value = values[i]
		sum_values += value
		if value < min_value
			min_value = value
			min_index = i
		elseif value > max_value
			max_value = value
			max_index = i
		end
	end
	mean_value = sum_values / length(values)
	mid_index = (1 + length(values)) ÷ 2
	# direct sorting of 3 indices without allocating arrays
	if min_index <= mid_index
		if mid_index <= max_index
			# min, mid, max order
			return SVector{3, Float64}(min_value, mean_value, max_value)
		elseif min_index <= max_index
			# min, max, mid order
			return SVector{3, Float64}(min_value, max_value, mean_value)
		else
			# max, min, mid order
			return SVector{3, Float64}(max_value, min_value, mean_value)
		end
	else
		if min_index <= max_index
			# mid, min, max order
			return SVector{3, Float64}(mean_value, min_value, max_value)
		elseif mid_index <= max_index
			# mid, max, min order
			return SVector{3, Float64}(mean_value, max_value, min_value)
		else
			# max, mid, min order
			return SVector{3, Float64}(max_value, mean_value, min_value)
		end
	end
end

WordStyle(::ESAX) = MultiWord{3}()



