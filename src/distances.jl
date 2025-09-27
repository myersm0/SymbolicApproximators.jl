
function Distances.evaluate(
		dist::Euclidean, 
		word1::Word{<:ContinuousApproximator, T}, 
		word2::Word{<:ContinuousApproximator, T}
	) where T
	n = length(word1.symbols)
	w = length(word2.symbols)
	n == w || throw(DimensionMismatch("Words must have same length"))
	# standard Euclidean distance on the PAA coefficients
	# todo: multiply by sqrt(n/w) for the compression rate factor (from paper eq. 4)
	return euclidean(word1.symbols, word2.symbols)
end

function Distances.evaluate(
		dist::Euclidean, word1::Word{SA, T}, word2::Word{SA, T}
	) where {SA <: SymbolicApproximator, T}
	n = length(word1)
	w = length(word2)
	n == w || throw(DimensionMismatch("Words must have same length"))
	keys1 = keys(word1)
	keys2 = keys(word2)
	approximator = word1.approximator[]
	total = 0.0
	for (i, j) in zip(keys1, keys2)
		d = symbol_distance(approximator, i, j)
		total += d^2
	end
	return sqrt(total)
end

function symbol_distance(approximator::SAX, i::Int, j::Int)
	α = alphabet(approximator)
	β = breakpoints(approximator)
	if abs(i - j) <= 1
       return 0.0
   end
	if i > j
		i, j = j, i
	end
	i += 1
	if i > length(β)
		i = length(β)
	end
	if j > length(β)
		j = length(β)
	end
   return abs(β[i] - β[j])
end


