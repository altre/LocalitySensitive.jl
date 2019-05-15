import DataStructures
import SimpleTraits
import Base.push!

struct MinHash
    salts :: Vector{UInt}
    """
    Construct minhash hash functions. 
    # Arguments
    - `threshold` the jaccard similarity above which documents are probably returned to be similar.
    """
    function MinHash(max_n_hashes=200)
        new([rand(UInt) for _ in 1:max_n_hashes])
    end
end

SimpleTraits.@traitfn function fingerprint(mh:: MinHash, shingles::::SimpleTraits.BaseTraits.IsIterator) 
    [minimum((hash(s, salt) for s in shingles)) for salt in mh.salts]
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
    julia> shingle = shinglerize(2)
    julia> collect(shingle("abcd"))
    ["ab","bc","cd"]
    ```
"""
function shingle(s::AbstractString; size = 3)::Vector{AbstractString}
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

mutable struct MinHashIndex
    threshold :: Float64
    tables :: Vector{DataStructures.DefaultDict{UInt, Vector{Int}}}
    bands :: Int
    rows :: Int
    current_index :: Int
    max_n_hashes :: Int
    function MinHashIndex(;minhash = MinHash(), threshold=0.5) 
        max_n_hashes  = length(minhash.salts)
        bands, rows = _get_partition(threshold, max_n_hashes)
        n_tables = bands * rows < max_n_hashes ? bands + 1 : bands
        tables = [DataStructures.DefaultDict{UInt, Vector{Int}}(Vector{Int}) for _ in 1:n_tables]
        new(threshold, tables, bands, rows, 1, max_n_hashes)
    end
end

function _get_partition(threshold, max_n_hashes)
    min_error = Inf
    opt = (0, 0)
    for b in 1 : max_n_hashes
        max_r = Int(floor(max_n_hashes / b))
        for r in 1 : max_r
            fp = _integrate(s -> 1 - (1 - s^r)^b, 0.0, threshold)
            fn = _integrate(s -> (1 - s^r)^b, threshold, 1.0)
            error = fp + fn
            if error < min_error
                min_error = error
                opt = (b, r)
            end
        end
    end
    return opt
end

function _integrate(f, a, b)
    p = 0.001
    area = 0.0
    x = a
    while x < b
        area = area + f(x+0.5*p)*p
        x = x + p
    end
    return area
end

"""
Push new signature to index.
"""
function push!(mhind:: MinHashIndex, signature::Vector{UInt})
    index = mhind.current_index
    hashes = _bands(mhind, signature)
    for (signature_part, table) in zip(hashes, mhind.tables)
        push!(table[signature_part], index)
    end
    mhind.current_index += 1;
end

function _bands(mhind, signature)
    (hash(signature[(i - 1) * mhind.rows + 1 : min(i * mhind.rows, mhind.max_n_hashes)]) for i in 1:length(mhind.tables))
end

"""
Find probably similar indices in index.
"""
function find_similar(mhind:: MinHashIndex, signature::Vector{UInt})::Vector{Int}
    indices = Vector{Vector{Int}}()
    hashes = _bands(mhind, signature)
    for (signature_part, table) in zip(hashes, mhind.tables)
        #TODO: remove double hash
        push!(indices, table[signature_part])
    end
    unique(vcat(indices...))
end

"""
Find all pairs of similar strings in index.
"""
function similar_pairs(mh:: MinHashIndex)::Vector{Tuple{Int, Int}}
    sims = Set{Tuple{Int, Int}}()
    for table in mh.tables
        for similars in values(table)
            for i in 1:length(similars), j in i + 1 : length(similars)
                push!(sims, (similars[i],similars[j]))
            end
        end
    end
    collect(sims)
end
