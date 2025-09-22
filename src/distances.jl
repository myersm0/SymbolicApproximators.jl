
function distance(d::SAX, symbols1::Vector{Char}, symbols2::Vector{Char})
	length(symbols1) == length(symbols2) || error("Symbol sequences must have same length")
	if !isnothing(d.distance_table)
		# use precomputed lookup table
		dist_sum = 0.0
		for i in 1:length(symbols1)
			idx1 = Int(symbols1[i] - 'a') + 1
			idx2 = Int(symbols2[i] - 'a') + 1
			dist_sum += d.distance_table[idx1, idx2]^2
		end
		return sqrt(dist_sum)
	else
		# fall back to computation
		dist_sum = 0.0
		for i in 1:length(symbols1)
			dist_sum += distance(symbols1[i], symbols2[i], d.breakpoints)^2
		end
		return sqrt(dist_sum)
	end
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

function distance(d::OrdinalDiscretizer, symbols1::Vector{Char}, symbols2::Vector{Char})
	# simple Hamming distance for ordinal
	length(symbols1) == length(symbols2) || error("Symbol sequences must have same length")
	return sum(symbols1 .!= symbols2)
end

function distance(d::DeltaDiscretizer, symbols1::Vector{Char}, symbols2::Vector{Char})
	# Hamming distance for delta (could be refined)
	length(symbols1) == length(symbols2) || error("Symbol sequences must have same length")
	return sum(symbols1 .!= symbols2)
end


