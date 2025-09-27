using SymbolicApproximators
using StatsBase
using Random
using Test

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
		# paper claims normalized data produces equiprobable symbols
		Random.seed!(42)
		n_samples = 10000
		data = sort(randn(n_samples))  # Standard normal data
		sax = SAX(1000, 5)  # 5 symbols
		symbols = encode(sax, data)
		counts = countmap(symbols)
		# each symbol should appear ~20% of the time (1000/5 = 200)
		expected_count = length(symbols) / 5
		for (sym, count) in counts
			# tolerance for randomness
			@test abs(count - expected_count) < expected_count * 0.1
		end
	end

end
