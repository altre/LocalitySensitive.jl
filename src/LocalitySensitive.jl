module LocalitySensitive
    include("minhash.jl")
    
    export MinHash, push!, find_similar, similar_pairs, fingerprint, estimate_jaccard, MinHashIndex, shingle, fingerprint_all
end # module
