# SymbolicDiscretizers
A Julia package to transform continuous signals (such as time series) into symbolic representations, while preserving meaningful patterns and distances.

** This package is still under development. Check back soon. **

In particular we implement the following algorithms:
- Symbolic Aggregate approXimation (SAX) ([Lin et al 2003](https://www.cs.ucr.edu/~eamonn/SAX.pdf))
- Extended SAX (ESAX) (Lkhagva et al 2006)
- Indexable SAX (iSAX) ([Shieh et al 2008](https://www.cs.ucr.edu/~eamonn/iSAX.pdf))
- Trend Feature SAX (TFSAX) ([Yu et al 2019](https://arxiv.org/abs/1905.00421)
- ordinal patterns ([Bandt & Pompe 2002](https://pubmed.ncbi.nlm.nih.gov/12005759/))
- delta encodings (todo: find reference)

## Usage

```julia
using SymbolicDiscretizers

signal = sin.(range(0, 4π, length=100))

sax = SAX(10, 5)  # discretize using SAX with 10 segments, 5 symbols
symbols = discretize(sax, signal)
# Result: ['c', 'd', 'e', 'd', 'c', 'b', 'a', 'b', 'c', 'd']

# compute distance between two symbolic sequences
symbols2 = discretize(sax, cos.(range(0, 4π, length=100)))
dist = distance(sax, symbols, symbols2)
```

## Discretizers

### SAX
- Divides signal into segments (PAA)
- Maps to symbols using Gaussian breakpoints
- Guarantees lower-bounding distances
- Best for general time series, normalized patterns

### OrdinalDiscretizer
- Preserves rank ordering
- Two methods: `:rank` and `:quantile`
- Best when relative position matters more than absolute values

### DeltaDiscretizer
- Encodes differences between consecutive values
- Configurable order (1st, 2nd differences)
- Can chain with other discretizers
- Best for interval patterns or rate of change analysis

## Parameter Selection
- **`nsegments`** (SAX): Start with length(data)/10 to length(data)/20
- **`alphabet_size`**: 
  - 3-4: Coarse patterns
  - 5-8: Balance
  - 10+: Finer detail, less compression
- **`order`** (Delta): 1 for velocity, 2 for acceleration

## References
- Lin, J., Keogh, E., Lonardi, S., & Chiu, B. (2003). A symbolic representation of time series, with implications for streaming algorithms. DMKD.

## License

MIT


[![Build Status](https://github.com/myersm0/SymbolicDiscretizers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/myersm0/SymbolicDiscretizers.jl/actions/workflows/CI.yml?query=branch%3Amain)
