using SymbolicApproximators
using StatsBase
using Distances
using GLMakie

signal1 = (sin.(range(0, 4π, length = 100)) |> x -> (x .- mean(x)) ./ std(x))
signal2 = (cos.(range(0, 4π, length = 100)) |> x -> (x .- mean(x)) ./ std(x))
true_distance = euclidean(signal1, signal2)

parameters = [(5, 5), (10, 10), (25, 25), (50, 50), (100, 100)]
results = map(parameters) do (w, α)
	sax = SAX(w, α)
	sax_distance = evaluate(MinDist(), encode(sax, signal1), encode(sax, signal2))
	bits = w * log2(α)  # total bits needed to store the SAX representation
	error = 100 * (true_distance - sax_distance) / true_distance
	(bits = bits, error = error)
end

x = [r.bits for r in results]
y = [r.error for r in results]

fig = Figure()
ax = Axis(
	fig[1,1], 
	xlabel = "Representation size (bits)", ylabel = "Approximation error (%)", 
	title = "SAX distance convergence", 
	xscale = log10
)
lines!(ax, x, y)
scatter!(ax, x, y, markersize = 15)
ylims!(ax, (0, 50))




