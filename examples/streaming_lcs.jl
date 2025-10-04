
# This demo streams random signal data, converts sliding windows (200 timepoints in length)
# into lower-dimensional symbolic SAX representations (word size 20), and continually tracks 
# and plots the longest repeated substring found so far.

# As each new SAX word arrives, we process symbols one-by-one: if a symbol continues a pattern 
# already seen in history, we extend the current match; otherwise, we saves the match to history 
# (updating the best-found length so far if needed) and start fresh.

# The three figure panels will show, from top to bottom: 
# 1. the raw incoming signal
# 2. the current SAX word (symbolic approximation of the signal)
# 3. the longest repeated pattern discovered across the entire stream so far


using SymbolicApproximators
using StatsBase
using DataStructures
using Observables
using GLMakie
using SuffixAutomata

w = 20
a = 10
window_len = 200

fps = 30
Δt = 1 / fps

sax = SAX(w, a)
buf = CircularBuffer{Float64}(window_len)
dest = Vector{Int}(undef, w)

signal = Observable(Vector{Float64}(undef, window_len))
current_word = Observable(Vector{Int}(undef, w))
lcs_match = Observable(Int[])
lcs_length = Observable(0)

fig = Figure(; size = (800, 600))

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

history = SuffixAutomaton(Int[])
current_match = Int[]
best_match = Int[]

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
				if length(history) > 0 && occursin(candidate, history)
					push!(current_match, symbol)
				else
					if length(current_match) > length(best_match)
						global best_match = copy(current_match)
						lcs_match[] = best_match
						lcs_length[] = length(best_match)
					end
					append!(history, current_match)
					empty!(current_match)
					if occursin([symbol], history)
						push!(current_match, symbol)
					else
						push!(history, symbol)
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


`
