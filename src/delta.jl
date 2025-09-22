
struct DeltaDiscretizer <: SymbolicDiscretizer
	alphabet_size::Int
	order::Int
	base_discretizer::Union{Nothing, SymbolicDiscretizer}
	function DeltaDiscretizer(alphabet_size::Int; order::Int=1, base::Union{Nothing, SymbolicDiscretizer}=nothing)
		order > 0 || error("order must be positive")
		new(alphabet_size, order, base)
	end
end

function discretize(d::DeltaDiscretizer, series::Vector{Float64})
	diff_series = series
	for _ in 1:d.order
		diff_series = diff(diff_series)
	end
	isnothing(d.base_discretizer) || return discretize(d.base_discretizer, diff_series)
	isempty(diff_series) && return Char[]
	min_val, max_val = extrema(diff_series)
	range_val = max_val - min_val
	if range_val < 1e-10
		middle_symbol = Char('a' + div(d.alphabet_size, 2))
		return fill(middle_symbol, length(diff_series))
	end
	symbols = Vector{Char}(undef, length(diff_series))
	for i in 1:length(diff_series)
		normalized = (diff_series[i] - min_val) / range_val
		bin_idx = min(floor(Int, normalized * d.alphabet_size), d.alphabet_size - 1)
		symbols[i] = Char('a' + bin_idx)
	end
	return symbols
end

function distance(d::DeltaDiscretizer, symbols1::Vector{Char}, symbols2::Vector{Char})
	# Hamming distance for delta (could be refined)
	length(symbols1) == length(symbols2) || error("Symbol sequences must have same length")
	return sum(symbols1 .!= symbols2)
end


