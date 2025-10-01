
struct MinDist <: Metric end
const MINDIST = MinDist  # spell it the way the papers do, if desired

"""
    mindist(word1, word2)

Evaluate the MINDIST between two words.
"""
mindist(word1, word2) = evaluate(MinDist(), word1, word2)

# todo: factor out checking of n, word size, etc

function Distances.evaluate(
		dist::MinDist, 
		word1::Word{<:ContinuousApproximator}, 
		word2::Word{<:ContinuousApproximator}
	)
	word_size(word1) == word_size(word2) || error("Words must have same length")
	word1.n == word2.n || error("Words must have the same n")
	alphabet(word1) == alphabet(word2) || error("Words must have the same alphabet")
	return euclidean(values(word1), values(word2)) * sqrt(compression_rate(word1))
end

function Distances.evaluate(
		dist::MinDist, word1::Word{SA, T, 1}, word2::Word{SA, T, 1}
	) where {SA <: SymbolicApproximator, T}
	word_size(word1) == word_size(word2) || error("Words must have same length")
	word1.n == word2.n || error("Words must have same n")
	alphabet(word1) == alphabet(word2) || error("Words must have the same alphabet")
	keys1 = keys(word1)
	keys2 = keys(word2)
	approximator = word1.approximator[]
	total = 0.0
	for (i, j) in zip(keys1, keys2)
		d = symbol_distance(approximator, i, j)
		total += d^2
	end
	return sqrt(total) * sqrt(compression_rate(word1))
end

function Distances.evaluate(
		dist::MinDist, word1::Word{SA, T, W}, word2::Word{SA, T, 1}
	) where {SA <: SymbolicApproximator, T, W}
	error("Distance function not yet implemented for MultiWords")
end

function symbol_distance(approximator::SAX, i::Int, j::Int)
	β = breakpoints(approximator)
	if abs(i - j) <= 1
       return 0.0
   end
	if i > j
		i, j = j, i
	end
	lo = 0
	hi = length(β)
	i = clamp(i, lo:hi)
	j = clamp(j, lo:hi)
   return β[j] - β[i+1]
end


