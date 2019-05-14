module LocalitySensitive
    import Base.push!
    import DataStructures
    import StatsBase
    import SimpleTraits

    include("minhash.jl")
    
    export MinHash, push!, find_similar, similar_pairs, fingerprint, shinglerize, estimate_jaccard, MinHashIndex
end # module
