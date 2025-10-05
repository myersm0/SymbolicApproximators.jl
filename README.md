# SymbolicApproximators

## Overview
A Julia package implementing **Symbolic Aggregate approXimation (SAX)** and related methods for symbolic discretization and dimension reduction. These techniques allow you to take continuous-valued sequences (including but not necessarily limited to time series) and then use things like string algorithms and machine learning algorithms for categorical data.

The following algorithms are implemented or under development:
| |Algorithm|
|-|:-------------------------------------------------------------------------------------------------|
|☑|Piecewise Aggregate Approximation (PAA) ([Keogh et al 2000](https://www.cs.ucr.edu/~eamonn/SAX.pdf) and [Yi et al 2000](https://dl.acm.org/doi/10.5555/645926.671689))|
|☑|Symbolic Aggregate approXimation (SAX) ([Lin et al 2003](https://www.cs.ucr.edu/~eamonn/SAX.pdf))|
|☑|Extended SAX (ESAX) (Lkhagva et al 2006)|
|☐|Indexable SAX (iSAX) ([Shieh et al 2008](https://www.cs.ucr.edu/~eamonn/iSAX.pdf))|
|☐|Trend Feature SAX (TFSAX) ([Yu et al 2019](https://arxiv.org/abs/1905.00421)|
|☐|ordinal patterns ([Bandt & Pompe 2002](https://pubmed.ncbi.nlm.nih.gov/12005759/))|
|☐|delta encodings (todo: find reference)|

Note that the first one, PAA, is actually a continuous rather than a symbolic representation, but we include it in the package because all the SAX variants depend on it.

Coming soon there will be additional functionality such as:
- other functions depending on algorithm, such as numerosity reduction and permutation entropy

## Installation
```julia
using Pkg
Pkg.add("SymbolicApproximators")
```

## Usage

### Basic workflow
Usage revolves around this basic workflow:
1. **Preprocess your data**. See preprocessing section below. This package does not currently provide preprocessing utilities like z-normalization, but most algorithms assume it has been done.
2. Define a **`SymbolicApproximator`** with integer arguments for:
    - **word size**, i.e. the number of segments to use, or in other words the desired output length
    - **alphabet size**, also called cardinality
3. Pass that approximator and your normalized data into the function **`encode(::SymbolicApproximator, ::AbstractVector)`**.
4. Your output will be a **`Word`** composed of instances of the symbol set defined in the approximator.

```julia
julia> using SymbolicApproximators
julia> using StatsBase

julia> signal = sin.(range(0, 4π, length=100))
julia> normalized = (signal .- mean(signal)) ./ std(signal, corrected = false)

julia> approximator = SAX(10, 5)  # discretize using SAX with 10 segments, 5 symbols
julia> word = encode(approximator, normalized)
SAX Word: "decabdecab"
```

### An important note about preprocessing
Some algorithms expect preprocessed inputs. Specifically, SAX and variants expect the data to be normalized with a mean of 0, standard deviation of 1. It's left to the user to handle such preprocessing where necessary, mainly because there are a number of ways you can do it, depending on your situation: maybe your data already happens to be normally distributed in this manner, or maybe you have streaming data and need to do online normalization, etc.

It's also left to the user whether and how the time series is to be divided into "subsequences" (i.e. a scheme where you will have multiple _words_ per time series), and, if so, whether to normalize each subsequence independently. The authors of SAX do so: "we normalize each time series (including subsequences) to have a mean of zero and a standard deviation of one." In addition they advise: "if the subsequence contains only one value, the standard deviation is not defined. More troublesome is the case where the subsequence is _almost_ constant [...]. We can easily deal with this problem, if the standard deviation of the sequence before normalization is below an epsilon ϵ, we simply assign the entire word to the middle-ranged alphabet."

### Extracting content from a `Word`
The result of an `encode()` or `approximate()` call will be a `Word` struct. Internally, a `Word` stores integer indices into the alphabet of symbols, rather than storing the symbols themselves. Two accessors are provided to retrieve the contents:
- use **`keys(w::Word)`** to get the _integer indices_
- use **`values(w::Word)`** to get a vector of the _symbols_
    - note that the symbol values are generated from integers upon demand when you call `values()`, so, if execution time is _very_ important to you, you may find it sufficient to just use `keys()` instead

```julia
julia> keys(word)
10-element Vector{Int64}:
 4
 5
 ⋮
 2

julia> values(word)
10-element Vector{Char}:
 'd': ASCII/Unicode U+0064 (category Ll: Letter, lowercase)
 'e': ASCII/Unicode U+0065 (category Ll: Letter, lowercase)
 ⋮
 'b': ASCII/Unicode U+0062 (category Ll: Letter, lowercase)

julia> width(word)
1
```

Note the last operation there, `width`. Most algorithms produce `Word`s with a single symbol per position, but some (such as ESAX) represent multiple symbols in each position. We'll refer to this as a word's _width_. ESAX, for example, has a width of 3:
```julia
julia> esax = ESAX(10, 5)
ESAX(10, 5)

julia> word = encode(esax, normalized)
ESAX Word: SVector{3, Char}[['c', 'd', 'e'], ['e', 'e', 'e'], …, ['a', 'b', 'c']]

julia> values(word)
10-element Vector{SVector{3, Char}}:
 ['c', 'd', 'e']
 ['e', 'e', 'e']
 ⋮
 ['a', 'b', 'c']

julia> width(word)
3
```

### Distances

The [Distances.jl](https://github.com/JuliaStats/Distances.jl) framework is extended here to implement the MINDIST Euclidean-like distance algorithm for SAX and PAA words:

``` julia
signal1 = (sin.(range(0, 4π, length=100)) |> x -> (x .- mean(x)) ./ std(x))
signal2 = (cos.(range(0, 4π, length=100)) |> x -> (x .- mean(x)) ./ std(x))

sax = SAX(50, 10)
word1 = encode(sax, signal1)
word2 = encode(sax, signal2)

# Euclidean-like "MINDIST" between the approximated time series:
evaluate(MinDist(), word1, word2) # 11.4

# or, equivalently:
mindist(word1, word2)             # 11.4

# true Euclidean distance between the original time series:
euclidean(signal1, signal2)       # 14.1

# if we use SAX with a much higher alphabet size, we can get very close
# to the true Euclidean distance:
sax = SAX(50, 250)
word1 = encode(sax, signal1)
word2 = encode(sax, signal2)
mindist(word1, word2)             # 13.9

```

The MINDIST between two SAX or PAA words lower-bounds the true Euclidian distance, and will closely approach that distance (at reduced representational cost) as we increase
the paramaters for word size and/or alphabet size. The example file [convergence.jl](https://github.com/myersm0/SymbolicApproximators.jl/blob/main/examples/convergence.jl) in this repo demonstates that.

### Advanced usage for streaming
For streaming or real-time processing applications, you can use the new mutating `encode!()` function to avoid allocations by reusing a pre-allocated destination array:
```julia
using SymbolicApproximators
using StatsBase
import SymbolicApproximators: windows

signal = randn(1000)
subsequences = windows(signal, 100)  # generate sliding windows of size 100

# define a preprocessing function; in practice this might be more complex
znormalize(x) = (x .- mean(x)) ./ std(x)

sax = SAX(10, 5)
n_words = length(subsequences)
dest = Matrix{Int}(undef, 10, n_words)
words = [
    encode!(sax, view(dest, :, i), znormalize(subseq))
    for (i, subseq) in enumerate(subsequences)
]
```

The contents of each word above will be a _view_ into the respective column of the destination array. This pattern may reduce allocations and improve memory contiguity. Depending on your use case, you may not even need the `Word` structs themselves, and the `dest` matrix may have all you need. (Specifically, it will contain the (_zero-based_) integer indices into the symbolic alphabet.)

See [streaming_basic.jl](https://github.com/myersm0/SymbolicApproximators.jl/blob/main/examples/streaming_basic.jl) for a complete example with real-time visualization, and [streaming_lcs.jl](https://github.com/myersm0/SymbolicApproximators.jl/blob/main/examples/streaming_lcs.jl) for finding repeated patterns in streaming SAX words.

## License

MIT


[![Build Status](https://github.com/myersm0/SymbolicApproximators.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/SymbolicApproximators.jl/actions/workflows/CI.yml?query=branch%3Amain)
