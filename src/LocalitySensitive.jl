module LocalitySensitive
    include("minhash.jl")
    include("simhash.jl")
    
    export MinHash, push!, find_similar, similar_pairs, fingerprint, estimate_jaccard, LSHIndex, shingle, fingerprint_all
    export SimHash, SimHashIndex, estimate_cosine
end # module
