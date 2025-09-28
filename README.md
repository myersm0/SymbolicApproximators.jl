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
- rolling window/segment functions to enable operations on streaming data
- other functions depending on algorithm, such as numerosity reduction and permutation entropy

## Usage

### Basic workflow
Usage revolves around this basic workflow:
1. Define a **`SymbolicApproximator`** with integer arguments for:
    - **word size**, i.e. the number of segments to use, or in other words the desired output length or dimensionality.
    - **alphabet size**, also called cardinality. This refers to the number of symbols (`Char`s by default) that you want to use.
        - Alternatively instead of specifying an alphabet _length_, you can directly pass in the _symbol set_ that you want to use in your output word. For example, `SAX(10, -2:2)` will give you an alphabet size of 5 where the symbols will be the integers -2 through 2 inclusive.
2. Pass that approximator and your data (presumably normalized -- see below) into the function **`encode(::SymbolicApproximator, ::AbstractVector)`**. Or if you prefer, use `approximate()` which is an alias for `encode()`.
3. Your output will be a **`Word`** composed of instances of the symbol set defined in the approximator.
    - A `Word` holds a representation of your encoded output, along with a reference to the `SymbolicApproximator` that generated it (to inform things like the computation of distance between two `Word`s, for example, which may be algorithm-dependent).


Note that some algorithms expect preprocessed inputs. Specifically, SAX and variants expect the data to be normalized with a mean of 0, standard deviation of 1. It's left to the user to handle such preprocessing where necessary, mainly because there are a number of ways you can do it, depending on your situation: maybe your data already happens to be normally distributed in this manner, or maybe you have streaming data and need to do online normalization, etc.

```julia
julia> using SymbolicApproximators
julia> using StatsBase

julia> signal = sin.(range(0, 4π, length=100))
julia> normalized = (signal .- mean(signal)) ./ std(signal, corrected = false)

julia> approximator = SAX(10, 5)  # discretize using SAX with 10 segments, 5 symbols
julia> word = encode(approximator, normalized)
SAX Word: "decabdecab"
```

### Extracting content from a `Word`
Internally, a `Word` stores integer indices into the alphabet of symbols, rather than storing the symbols themselves. Two accessors are provided to retrieve the contents:
- use **`keys(w::Word)`** to get the _integer indices_
- use **`values(w::Word)`** to get a vector of the _symbols_
    - note that the symbol values are mapped from integers upon demand when you call `values()`, so, if execution time is very important to you, you may find it sufficient to just use `keys()` instead

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

Note the last operation there, `width`. Most algorithms produce `Word`s with a single symbol per position, but some (such as ESAX) represent multiple symbols in each position. We'll refer to this as a word's _width_. ESAX, for example, has a `width` of 3:
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
julia> signal1 = (sin.(range(0, 4π, length=100)) |> x -> (x .- mean(x)) ./ std(x))
julia> signal2 = (cos.(range(0, 4π, length=100)) |> x -> (x .- mean(x)) ./ std(x))
julia> sax = SAX(50, 10)
julia> word1 = encode(sax, signal1)
julia> word2 = encode(sax, signal2)

julia> evaluate(MinDist(), word1, word2)
11.40858928469606

# or, equivalently:
julia> mindist(word1, word2)
11.40858928469606

julia> euclidean(signal1, signal2)
14.071247279470288
```

The MINDIST between two SAX or PAA words lower-bounds the true Euclidian distance, and will very closely approach that distance (at reduced representational cost) as we increase
the paramaters for word size and/or alphabet size. The example file [convergence.jl](https://github.com/myersm0/SymbolicApproximators.jl/blob/main/examples/convergence.jl) in this repo demonstates that.

## License

MIT


[![Build Status](https://github.com/myersm0/SymbolicApproximators.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/SymbolicApproximators.jl/actions/workflows/CI.yml?query=branch%3Amain)
