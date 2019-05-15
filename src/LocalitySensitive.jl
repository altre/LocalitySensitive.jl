module LocalitySensitive
    include("minhash.jl")
    
    export MinHash, push!, find_similar, similar_pairs, fingerprint, shinglerize, estimate_jaccard, MinHashIndex
end # module
