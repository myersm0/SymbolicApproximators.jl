# SymbolicApproximators
A Julia package implementing Symbolic Aggregate approXimation (SAX) and related methods for symbolic time series discretization and dimension reduction.

In particular we have implemented (or are in the process of doing so) the following algorithms:
| |Algorithm|
|-|:------------------------------------------------------------------------------------------------|
|☑|Piecewise Aggregate Approximation (PAA) (todo: find reference)|
|☑|Symbolic Aggregate approXimation (SAX) ([Lin et al 2003](https://www.cs.ucr.edu/~eamonn/SAX.pdf))|
|☐|Extended SAX (ESAX) (Lkhagva et al 2006)|
|☐|Indexable SAX (iSAX) ([Shieh et al 2008](https://www.cs.ucr.edu/~eamonn/iSAX.pdf))|
|☐|Trend Feature SAX (TFSAX) ([Yu et al 2019](https://arxiv.org/abs/1905.00421)|
|☐|ordinal patterns ([Bandt & Pompe 2002](https://pubmed.ncbi.nlm.nih.gov/12005759/))|
|☐|delta encodings (todo: find reference)|

Note that the first one, PAA, is actually a continuous rather than a symbolic representation, but we include it in the package because all the SAX variants depend on it.

## Usage
Note that some algorithms expect preprocessed inputs. Specifically, SAX and variants expect the data to be normalized with a mean of 0, standard deviation of 1. We deliberately leave it to the user to handle such preprocessing where necessary, mainly because there are a number of ways you can do it: maybe your data already happens to be normally distributed in this manner, or maybe you have streaming data and need to do online normalization, etc.

```julia
using SymbolicApproximators
using StatsBase

signal = sin.(range(0, 4π, length=100))
normalized = (signal .- mean(signal)) ./ std(signal, corrected = false)

approximator = SAX(10, 5)  # discretize using SAX with 10 segments, 5 symbols
symbols = encode(approximator, normalized)
# Result: ['d', 'e', 'c', 'a', 'b', 'd', 'e', 'c', 'a', 'b']
```

Operations revolve around this basic workflow:
1. Define a `SymbolicApproximator` with integer arguments for _word size_ (i.e., the number of segments) and _alphabet size_ (also called cardinality), respectively. For example, `SAX(10, 5)`.
  - Alternatively instead of specifying an alphabet _length_, you can directly pass in the _symbol set_ that you want to use in your output word. For example, `SAX(10, -2:2)` or equivalently `SAX(10, [-2, 1, 0, 1, 2])` will give you an alphabet size of 5 where the symbols will be the numbers -2 through 2 inclusive.)
2. Pass that approximator and your (presumably normalized) data into the `encode()` function.
3. Your output will be a `Word <: AbstractVector` composed of instances of the symbol set defined in the approximator.

Coming soon there will be additional functionality such as:
- more algorithms implemented
- distance functions (using [Distances.jl](https://github.com/JuliaStats/Distances.jl))
- other functions depending on algorithm, such numerosity reduction and permutation entropy

## License

MIT


[![Build Status](https://github.com/myersm0/SymbolicApproximators.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/SymbolicApproximators.jl/actions/workflows/CI.yml?query=branch%3Amain)
