
function numerosity_reduction(symbols::Vector{Char})
	# run-length encoding for consecutive identical symbols
	vals, lens = rle(symbols)
	return vals, lens
end

function expand_numerosity(vals::Vector{Char}, lens::Vector{Int})
	return inverse_rle(vals, lens)
end

function reconstruct(sax::SAX, symbols::Vector{Char}, original_length::Int)
	paa_values = zeros(length(symbols))
	
	for i in 1:length(symbols)
		idx = Int(symbols[i] - 'a')
		if idx == 0
			paa_values[i] = sax.breakpoints[1] - 0.5
		elseif idx >= length(sax.breakpoints)
			paa_values[i] = sax.breakpoints[end] + 0.5
		else
			paa_values[i] = (sax.breakpoints[idx] + sax.breakpoints[idx+1]) / 2
		end
	end
	
	segment_length = original_length / length(symbols)
	reconstructed = Float64[]
	
	for i in 1:length(symbols)
		repeat_count = round(Int, segment_length)
		if i == length(symbols)
			repeat_count = original_length - length(reconstructed)
		end
		append!(reconstructed, fill(paa_values[i], repeat_count))
	end
	
	return reconstructed[1:original_length]
end

