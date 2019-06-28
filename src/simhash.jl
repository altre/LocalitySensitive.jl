import SparseArrays
import BKTrees
import Base.push!
import SparseArrays

struct SimHash
    salt::UInt
    SimHash(seed = rand(UInt)) = new(seed)
end

function fingerprint(sh::SimHash, features_weights)
    fingerprints = zeros(64)
    for (feature, weight) in features_weights
        hash_ = hash(feature, sh.salt)
        @inbounds @simd for i in 1:64
            if (hash_ & _bitmask(i)) != 0
                fingerprints[i] += weight
            else
                fingerprints[i] -= weight
            end
        end
    end
    signature = UInt(0)
    for (i, finger) in enumerate(fingerprints)
        if finger > 0
            signature |= _bitmask(i)
        end
    end
    signature
end

function _bitmask(i)
    UInt(1) << (i - 1)
end

function fingerprint_all(sh::SimHash, data::SparseArrays.SparseMatrixCSC)
    n = size(data, 2)
    fingerprints = Vector{UInt64}(undef, n)
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

bk_distfunc(a,b) = hamming(a[2], b[2])

mutable struct SimHashIndex
    tree::BKTrees.BKTree
    current_index::Int
    distance_bits::Int
    SimHashIndex(cosine_similarity_threshold) = new(BKTrees.BKTree{Tuple{Int, UInt}}(bk_distfunc), 
            1, _cosine_to_hamming(cosine_similarity_threshold))
end

function _cosine_to_hamming(t)
    Int(round(acos(t)/π*64))
end

function push!(sh::SimHashIndex, fingerprint)
    BKTrees.add!(sh.tree, (sh.current_index, fingerprint))
    sh.current_index += 1;
end

function find_similar(sh::SimHashIndex, fingerprint)
    [t[2][1] for t in BKTrees.find(sh.tree, (0,fingerprint), sh.distance_bits; k=sh.current_index)]
end

function similar_pairs(shind)
    pairs = []
    nodes = [shind.tree.root]
    while !isempty(nodes)
        current = pop!(nodes)
        if current.item == nothing
            break
        end
        for c in values(current.children)
            push!(nodes, c)
        end
        similar = BKTrees.find(shind.tree, current.item, shind.distance_bits, k=shind.current_index)
        for (_, (index, _)) in similar
            if current.item[1] < index
                push!(pairs, (current.item[1], index))
            end
        end
    end
    pairs
end