# Examples

## convergence.jl
Shows how SAX approximation quality improves with increased representation size: the MINDIST between two SAX words converges to the true Euclidean distance as you increase the word size and alphabet size parameters.

## lcs.jl  
Finds longest common substring between SAX-encoded signals.

Encodes sine and cosine waves as SAX strings and uses a suffix automaton structure to find their longest common substring. The visualization shows both the original signals and their SAX representations, with the matching pattern highlighted.

## streaming_basic.jl
An example of real-time SAX encoding of streaming data.

A live visualization that continuously:
- Generates noisy sine wave samples
- Maintains a sliding window buffer  
- Encodes each window as a SAX word using pre-allocated memory
- Updates plots showing both raw signal and its symbolic representation

## streaming_lcs.jl
Online pattern discovery in streaming SAX words approximated from a random signal. This demo continuously searches for the longest repeated pattern across the entire history of observed SAX words.

It works incrementally - as each new symbol arrives, we either extend the current match (if that extended sequence exists in history) or terminate the current match and starts fresh.
