
# PAA (Piecewise Aggregate Approximation)
function paa(series::Vector{Float64}, nsegments::Int)
	n = length(series)
	segment_length = n / nsegments
	result = zeros(nsegments)
	for i in 1:nsegments
		start = round(Int, (i - 1) * segment_length + 1)
		stop = min(n, round(Int, i * segment_length))
		result[i] = mean(series[start:stop])
	end
	return result
end

function value_to_symbol(value::Float64, breakpoints::Vector{Float64})
	for (i, b) in enumerate(breakpoints)
		value < b && return Char('a' + i - 1)
	end
	return Char('a' + length(breakpoints))
end

function value_to_symbol_ordinal(value::Float64, breakpoints::Vector{Float64})
	for (i, b) in enumerate(breakpoints)
		value <= b && return Char('a' + i - 1)
	end
	return Char('a' + length(breakpoints))
end

function linspace_bins(start, stop, n)
	return collect(range(start, stop, length=n+1))[2:end]
end

function distance(s1::Char, s2::Char, breakpoints::Vector{Float64})
	s1 == s2 && return 0.0
	
	index1 = Int(s1 - 'a')
	index2 = Int(s2 - 'a')
	
	abs(index1 - index2) > 1 || return 0.0
	
	maxindex = max(index1, index2)
	minindex = min(index1, index2)
	
	if minindex == 0
		low_val = breakpoints[1]
	else
		low_val = breakpoints[minindex]
	end
	
	if maxindex > length(breakpoints)
		high_val = breakpoints[end]
	else
		high_val = breakpoints[maxindex]
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

