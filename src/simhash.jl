import SparseArrays
import Base.push!
import SparseArrays
import DataStructures

include("index_util.jl")

struct SimHash
    salt::UInt
    n_hashes :: Int
    SimHash(seed = rand(UInt), n_hashes=64) = new(seed, n_hashes)
end

function fingerprint(sh::SimHash, features_weights)
    fingerprints = zeros(sh.n_hashes)
    for (feature, weight) in features_weights
        hash_ = BitVector()
        hash_.chunks = [hash(feature, sh.salt) for _ in 1:ceil(sh.n_hashes / 64)]
        hash_.dims = (sh.n_hashes,)
        hash_.len = sh.n_hashes
        @inbounds @simd for i in hash_
            if i != 0
                fingerprints[i] += weight
            else
                fingerprints[i] -= weight
            end
        end
    end
    signature = BitVector(undef, sh.n_hashes)
    for (i, finger) in enumerate(fingerprints)
        signature[i] = finger > 0
    end
    signature
end

function _bitmask(i)
    UInt(1) << (i - 1)
end

function fingerprint_all(sh::SimHash, data::SparseArrays.SparseMatrixCSC)
    n = size(data, 2)
    fingerprints = Vector{BitVector}(undef, n)
    Threads.@threads for i in 1:n
        col = data[:,i]
        fingerprints[i] = fingerprint(sh, col)
    end
    fingerprints
end

function fingerprint(sh::SimHash, features_weights::SparseArrays.SparseVector)
    fingerprint(sh, zip(features_weights.nzind, features_weights.nzval))
end

"""
    Using the function P(h(a) = h(b)) = 1 - θ/π, return an estimate of the cosine distance.
"""
function estimate_cosine(a, b)
    cos((hamming(a,b) / 64) * π )
end

#TODO: Benchmark xor popcnt instead of distances packages
@inline function hamming(a, b)
    count_ones(xor(a,b))
end

mutable struct SimHashIndex <: LSHIndex
    tables :: Vector{DataStructures.DefaultDict{UInt64, Vector{Int}}}
    bands :: Int
    rows :: Int
    current_index :: Int
    masks :: Vector{UInt}
    function SimHashIndex(cosine_similarity_threshold)
        bands, rows = _get_partition(_cosine_to_probability(cosine_similarity_threshold), 64)
        tables = [DataStructures.DefaultDict{UInt64, Vector{Int}}(Vector{Int}) for _ in 1:bands]
        masks = [_mask((i - 1) * rows + 1 : min(i * rows, 64)) for i in 1:bands]
        new(tables, bands, rows, 1, masks)
    end
end

"""
Push new fingerprint to index.
"""
function push!(shind::SimHashIndex, fp::UInt)
    index = shind.current_index
    masks = shind.masks
    bands = fp .& masks
    for (fingerprint_part, table) in zip(bands, shind.tables)
        push!(table[fingerprint_part], index)
    end
    shind.current_index += 1;
end

"""
Find probably similar indices in index.
"""
function find_similar(shind::SimHashIndex, fp::UInt)::Vector{Int} 
    indices = Vector{Vector{Int}}()
    bands = fp .& shind.masks
    for (band, table) in zip(bands, shind.tables)
        push!(indices, table[band])
    end
    unique(vcat(indices...))
end

function _mask(r::UnitRange{Int})
    z = UInt(0)
    for i in r
        z |= (UInt(1) << (i - 1))
    end
    z
end

function _cosine_to_hamming(t)
    Int(round(acos(t)/π*64))
end

function _cosine_to_probability(t)
    1-acos(t)/π
end
