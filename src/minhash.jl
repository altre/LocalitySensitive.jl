import DataStructures
import SimpleTraits
import Random
import Base.push!

include("index_util.jl")

struct MinHash
    salts :: Vector{UInt}
    n_hashes :: Int
    """
    Construct minhash hash functions. 
    # Arguments
    - `threshold` the jaccard similarity above which documents are probably returned to be similar.
    """
    function MinHash(max_n_hashes=200, rng=Random.MersenneTwister())
        new([rand(rng, UInt) for _ in 1:max_n_hashes], max_n_hashes)
    end
end

SimpleTraits.@traitfn function fingerprint(mh:: MinHash, shngls::::SimpleTraits.BaseTraits.IsIterator) 
    [minimum((hash(s, salt) for s in shngls)) for salt in mh.salts]
end

function fingerprint_all(mh::MinHash, shingles)
    n_hashes = length(mh.salts)
    n_docs = length(shingles)
    all_fingerprints = Array{UInt, 2}(undef,n_hashes,n_docs) 
    Threads.@threads for i in 1:n_docs
        all_fingerprints[:,i] = fingerprint(mh, shingles[i])
    end
    [all_fingerprints[:,i] for i in 1:n_docs]
end

"""
    Compute MinHash fingerprint for string using the `shingle` function.
"""
function fingerprint(mh:: MinHash, str::AbstractString; shingle_size=3)::Vector{UInt}
    fingerprint(mh, shingle(str; size=shingle_size))
end

"""
    Cut a given string into shingles of size `size`.
    ```jldoctest
    julia> shingle("abcd"; size = 2)
    ["ab","bc","cd"]
    ```
"""
function shingle(s::AbstractString; size = 3)::Vector{String}
    if length(s) <= size
        return [s]
    end
    shingles = Vector{String}()
    for i in 1:(length(s) - size + 1)
        push!(shingles, s[thisind(s, i): thisind(s, i + size - 1)])
    end
    unique(shingles)
end

function estimate_jaccard(a::Vector{UInt}, b::Vector{UInt})
    sum(a .== b) / length(a)
end

mutable struct MinHashIndex <: LSHIndex
    threshold :: Float64
    tables :: Vector{DataStructures.DefaultDict{Vector{UInt}, Vector{Int}}}
    bands :: Int
    rows :: Int
    current_index :: Int
    max_n_hashes :: Int
    function MinHashIndex(;fingerprinter = MinHash(), threshold=0.5)
        max_n_hashes = fingerprinter.n_hashes
        bands, rows = _get_partition(threshold, max_n_hashes) 
        tables = [DataStructures.DefaultDict{Vector{UInt}, Vector{Int}}(Vector{Int}) for _ in 1:bands]
        new{Vector{UInt}}(threshold, tables, bands, rows, 1, max_n_hashes)
    end
end

"""
Push new fingerprint to index.
"""
function push!(mhind:: MinHashIndex, fingerprint::Vector{Int}) where {FingerprintType, HashVectorType}
    index = mhind.current_index
    hashes = _bands(mhind, fingerprint)
    for (fingerprint_part, table) in zip(hashes, mhind.tables)
        push!(table[fingerprint_part], index)
    end
    mhind.current_index += 1;
end

function _bands(mhind, fingerprint::Vector{Int})
    (fingerprint[(i - 1) * mhind.rows + 1 : min(i * mhind.rows, mhind.max_n_hashes)] for i in 1:length(mhind.tables))
end

"""
Find probably similar indices in index.
"""
function find_similar(mhind:: MinHashIndex{HashVectorType}, 
    fingerprint::HashVectorType)::Vector{Int} where {FingerprintType, HashVectorType}
    indices = Vector{Vector{Int}}()
    bands = _bands(mhind, fingerprint)
    for (band, table) in zip(bands, mhind.tables)
        push!(indices, table[band])
    end
    unique(vcat(indices...))
end

