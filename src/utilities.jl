
function compute_gaussian_breakpoints(alphabet_size)
	alphabet_size > 1 || error("alphabet_size must be at least 2")
	return quantile.(Normal(), (1:(alphabet_size-1)) ./ alphabet_size)
end

function build_distance_table(breakpoints::Vector{Float64})
	n = length(breakpoints) + 1
	table = zeros(n, n)
	for r in 1:n
		for c in 1:n
			if abs(r - c) <= 1
				table[r, c] = 0.0
			else
				max_idx = max(r, c)
				min_idx = min(r, c)
				if min_idx == 1
					low_val = breakpoints[1]
				else
					low_val = breakpoints[min_idx - 1]
				end
				if max_idx > length(breakpoints)
					high_val = breakpoints[end]
				else
					high_val = breakpoints[max_idx - 1]
				end
				table[r, c] = abs(high_val - low_val)
			end
		end
	end
	return table
end

# PAA (Piecewise Aggregate Approximation)
function paa(series::Vector{Float64}, nsegments::Int)
	n = length(series)
	segment_length = n / nsegments
	result = zeros(nsegments)
	for i in 1:nsegments
		start_idx = round(Int, (i - 1) * segment_length + 1)
		end_idx = round(Int, i * segment_length)
		end_idx = min(end_idx, n)
		result[i] = mean(series[start_idx:end_idx])
	end
	return result
end

function value_to_symbol(value::Float64, breakpoints::Vector{Float64})
	for i in 1:length(breakpoints)
		if value < breakpoints[i]
			return Char('a' + i - 1)
		end
	end
	return Char('a' + length(breakpoints))
end

function value_to_symbol_ordinal(value::Float64, breakpoints::Vector{Float64})
	for i in 1:length(breakpoints)
		if value <= breakpoints[i]
			return Char('a' + i - 1)
		end
	end
	return Char('a' + length(breakpoints))
end

function linspace_bins(start, stop, n)
	return collect(range(start, stop, length=n+1))[2:end]
end

function distance(s1::Char, s2::Char, breakpoints::Vector{Float64})
	s1 == s2 && return 0.0
	
	idx1 = Int(s1 - 'a')
	idx2 = Int(s2 - 'a')
	
	abs(idx1 - idx2) > 1 || return 0.0
	
	max_idx = max(idx1, idx2)
	min_idx = min(idx1, idx2)
	
	if min_idx == 0
		low_val = breakpoints[1]
	else
		low_val = breakpoints[min_idx]
	end
	
	if max_idx > length(breakpoints)
		high_val = breakpoints[end]
	else
		high_val = breakpoints[max_idx]
	end
	
	return abs(high_val - low_val)
end

function numerosity_reduction(symbols::Vector{Char})
	# run-length encoding for consecutive identical symbols
	vals, lens = rle(symbols)
	return vals, lens
end

function expand_numerosity(vals::Vector{Char}, lens::Vector{Int})
	return inverse_rle(vals, lens)
end

