
struct SAX <: SymbolicDiscretizer
	nsegments::Int
	alphabet_size::Int
	normalize::Bool
	breakpoints::Vector{Float64}
	distance_table::Union{Nothing, Matrix{Float64}}
	function SAX(nsegments::Integer, alphabet_size::Integer; normalize = true, use_table = true)
		breakpoints = compute_gaussian_breakpoints(alphabet_size)
		distance_table = use_table ? build_distance_table(breakpoints) : nothing
		new(nsegments, alphabet_size, normalize, breakpoints, distance_table)
	end
end

function discretize(d::SAX, series::Vector{Float64})
	n = length(series)
	# todo: make a `normalize()` fn
	if d.normalize
		μ = mean(series)
		σ = std(series, corrected=false)
		# handle constant series
		if σ < 1e-10
			middle_symbol = Char('a' + div(d.alphabet_size, 2))
			return fill(middle_symbol, d.nsegments)
		end
		normalized = (series .- μ) ./ σ
	else
		normalized = series
	end
	paa_values = paa(normalized, d.nsegments)
	symbols = [value_to_symbol(paa_values[i], d.breakpoints) for i in 1:d.nsegments]
	return symbols
end

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

