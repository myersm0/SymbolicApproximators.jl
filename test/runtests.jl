using SymbolicApproximators
using StatsBase
using Distances
using Distributions
using Random
using Test

znormalize(x) = (x .- mean(x)) ./ std(x)

@testset "SAX paper verification tests" begin
	
	@testset "Gaussian breakpoints (table 3 from paper)" begin
		# verify breakpoints match Table 3 from the SAX paper
		expected_breakpoints = Dict(
			3 => [-0.43, 0.43],
			4 => [-0.67, 0.0, 0.67],
			5 => [-0.84, -0.25, 0.25, 0.84],
			6 => [-0.97, -0.43, 0.0, 0.43, 0.97],
			7 => [-1.07, -0.57, -0.18, 0.18, 0.57, 1.07],
			8 => [-1.15, -0.67, -0.32, 0.0, 0.32, 0.67, 1.15],
		)
		for (siz, expected) in expected_breakpoints
			sax = SAX(10, siz)
			@test all(abs.(breakpoints(sax)[2:end] .- expected) .< 0.01)
		end
	end
	
	@testset "Equiprobability of symbols" begin
		n_samples = 10000
		sax = SAX(1000, 5)
		# make equal-size bins spanning the real line quantiles
		quantiles = range(0, 1; length = n_samples+2)[2:end-1]
		data = quantile.(Normal(), quantiles)  # inverse CDF of Gaussian
		symbols = encode(sax, data)
		counts = countmap(symbols)
		expected_count = length(symbols) / 5
		for (sym, count) in counts
			@test count == expected_count
		end
	end

end

znormalize(x) = (x .- mean(x)) ./ std(x)

@testset "Distances" begin
	ts1 = sin.(range(0, 4π, length = 100)) |> znormalize
	ts2 = cos.(range(0, 4π, length = 100)) |> znormalize
	euclidean_dist = evaluate(Euclidean(), ts1, ts2)
	
	alphabet_sizes = word_sizes = [5, 10, 25, 50, 100]

	@testset "PAA convergence with word size" begin
		paa_dists = Float64[]
		paa_errors = Float64[]
		for w in word_sizes
			paa = PAA(w)
			paa_word1 = encode(paa, ts1)
			paa_word2 = encode(paa, ts2)
			paa_dist = evaluate(MinDist(), paa_word1, paa_word2)
			push!(paa_dists, paa_dist)
			push!(paa_errors, euclidean_dist - paa_dist)  # gap to true distance
			# PAA distance should be lower bound
			@test paa_dist <= euclidean_dist
		end
		# Statistical test for convergence using correlation
		# Errors should decrease with word size
		correlation = corspearman(word_sizes, paa_errors)
		@test correlation < -0.75  # strong negative correlation expected
	end
	
	@testset "SAX convergence with alphabet size" begin
		w = 10
		sax_errors = Float64[]
		for a in alphabet_sizes
			sax = SAX(w, a)
			sax_word1 = encode(sax, ts1)
			sax_word2 = encode(sax, ts2)
			sax_dist = evaluate(MinDist(), sax_word1, sax_word2)
			error = euclidean_dist - sax_dist
			push!(sax_errors, error)
			@test sax_dist <= euclidean_dist
		end
		# test convergence using correlation:
		correlation = corspearman(alphabet_sizes, sax_errors)
		@test correlation < -0.75
		# for log-linear relationship (common in convergence):
		log_alphabet = log.(alphabet_sizes)
		log_correlation = cor(log_alphabet, sax_errors)
		@test abs(log_correlation) > abs(correlation) || abs(correlation) > 0.8
	end
	
	@testset "SAX convergence with word size" begin
		a = 10
		sax_errors = Float64[]
		for w in word_sizes
			sax = SAX(w, a)
			sax_word1 = encode(sax, ts1)
			sax_word2 = encode(sax, ts2)
			sax_dist = evaluate(MinDist(), sax_word1, sax_word2)
			error = euclidean_dist - sax_dist
			push!(sax_errors, error)
			@test sax_dist <= euclidean_dist
		end
		@test corspearman(word_sizes, sax_errors) < -0.75
	end
	
	# test interaction of both parameters
	@testset "SAX combined parameter effects" begin
		results = []
		for w in word_sizes, a in alphabet_sizes
			sax = SAX(w, a)
			sax_word1 = encode(sax, ts1)
			sax_word2 = encode(sax, ts2)
			sax_dist = evaluate(MinDist(), sax_word1, sax_word2)
			error = euclidean_dist - sax_dist
			push!(results, (w = w, a = a, error = error))
		end
		ws = [r.w for r in results]
		as = [r.a for r in results]
		errors = [r.error for r in results]
		combined = ws .* as
		@test corspearman(combined, errors) < -0.7
	end
end	



