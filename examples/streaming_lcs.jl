
using SymbolicApproximators
using StatsBase
using DataStructures
using Observables
using GLMakie
using SuffixAutomata

w = 20
a = 10
window_len = 200

fps = 10
Δt = 1 / fps

sax = SAX(w, a)
buf = CircularBuffer{Float64}(window_len)
dest = Vector{Int}(undef, w)

signal = Observable(Vector{Float64}(undef, window_len))
current_word = Observable(Vector{Int}(undef, w))
lcs_match = Observable(Int[])
lcs_length = Observable(0)

fig = Figure(resolution = (800, 600))

ax1 = Axis(fig[1, 1], title = "Incoming signal")
lines!(ax1, 1:window_len, signal)
limits!(ax1, 1, window_len, -3.27, 3.27)

ax2 = Axis(fig[2, 1], title = "Current SAX word")
limits!(ax2, 1, w, -1, a)
ax2.yticks = (0:(a - 1), ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"])
scatter!(ax2, 1:w, current_word, markersize = 10)

ax3 = Axis(fig[3, 1], title = "Longest repeated substring")
lcs_plot = scatter!(
	ax3, lift(x -> 1:length(x), lcs_match), lcs_match, 
	markersize = 10, color = :purple
)
ax3.yticks = (0:(a - 1), ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"])
limits!(ax3, 1, w * 10, -1, a)

Label(fig[3, 2], lift(x -> "Length: $x", lcs_length), tellheight = false)

running = Ref(true)

history = Int[]
history_automaton = SuffixAutomaton(Int[])
current_match = Int[]
best_match = Int[]

# this will take a minute to start running
@async begin
	t = 0.0
	while running[]
		new_val = randn()
		push!(buf, new_val)
		if length(buf) == window_len
			normed = (buf .- mean(buf)) ./ std(buf, corrected = false)
			word = encode!(sax, dest, normed)
			word_vec = collect(keys(word))
			for symbol in word_vec
				candidate = vcat(current_match, symbol)
				if length(history) > 0 && occursin(candidate, history_automaton)
					push!(current_match, symbol)
				else
					if length(current_match) > length(best_match)
						global best_match = copy(current_match)
						lcs_match[] = best_match
						lcs_length[] = length(best_match)
					end
					append!(history, current_match)
					global history_automaton = SuffixAutomaton(history)
					empty!(current_match)
					if occursin([symbol], history_automaton)
						push!(current_match, symbol)
					else
						push!(history, symbol)
						global history_automaton = SuffixAutomaton(history)
					end
				end
			end
			signal[] = copy(buf)
			current_word[] = word_vec
		end
		sleep(1 / fps)
		t += 0.1
	end
end


