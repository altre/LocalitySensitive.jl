# LocalitySensitive

[![Build Status](https://travis-ci.com/altre/LocalitySensitive.jl.svg?branch=master)](https://travis-ci.com/altre/LocalitySensitive.jl)
[![Coveralls](https://coveralls.io/repos/github/altre/LocalitySensitive.jl/badge.svg?branch=master)](https://coveralls.io/github/altre/LocalitySensitive.jl?branch=master)

Implementations of Locality Sensitive Hashing schemes.
## MinHash
Implementation of [Minwise Independent Hashing](https://en.wikipedia.org/wiki/MinHash).

Example usage:
```julia
using LocalitySensitive

documents = readlines("resources/benchmark_data.csv")
mh = MinHash()
fingerprints = fingerprint_all(mh, [shingle(d, size=4) for d in documents])
estimate_jaccard(fingerprints[1], fingerprints[2])
mhind = MinHashIndex(minhash=mh, threshold=0.9)
for f in fingerprints[:]
    push!(mhind, f)
end
pairs = similar_pairs(mhind)
find_similar(mhind, fingerprint(mh, shingle(documents[210], size=4)))
```