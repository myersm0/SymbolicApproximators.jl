using SymbolicApproximators
using SuffixAutomata
using StatsBase
using GLMakie

signal1 = (sin.(range(0, 4π, length=100)) |> x -> (x .- mean(x)) ./ std(x))
signal2 = (cos.(range(0, 4π, length=100)) |> x -> (x .- mean(x)) ./ std(x))

sax = SAX(50, 10)
word1 = encode(sax, signal1)
word2 = encode(sax, signal2)
string1 = join(values(word1))
string2 = join(values(word2))

automaton = SuffixAutomaton(string2)
common_pattern, position = lcs(string1, automaton)
match_range = position:position+length(common_pattern)-1

fig = Figure(; size = (1200, 800))

ax1 = Axis(fig[1,1], title = "Original signals", ylabel = "Value")
lines!(ax1, signal1, label = "sin", linewidth=3)
lines!(ax1, signal2, label = "cos", linewidth=3)

ax2 = Axis(
	fig[2,1], 
	title = "SAX words", 
	subtitle = "(longest common substring: '$(common_pattern)' at position $position)", 
	xlabel = "Position", 
	ylabel = "Symbol"
)
scatter!(
	ax2, 1:50, [Int(c) - Int('a') for c in string1], 
	label="sin SAX", markersize = 8
)
scatter!(
	ax2, 1:50, [Int(c) - Int('a') for c in string2], 
	label="cos SAX", markersize = 8, marker = :diamond
)

vspan!(
	ax2, position, position+length(common_pattern)-1, 
	alpha = 0.2, color = :green
)
lines!(
	ax2, match_range, [Int(c) - Int('a') for c in common_pattern], 
	linewidth = 15, color = RGBA(0, 0, 0, 0.3), label = "LCS"
)

ax2.xticks = [1, position, 50]
axislegend(ax1)
axislegend(ax2)

save("sax_demo_lcs.png", fig)



