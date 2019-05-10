module LocalitySensitive
    import Base.push!
    import DataStructures

    struct MinHash
        threshold :: Float64
        documents :: Vector{AbstractString}
        shingle_size :: Int64
        n_bands :: Int64
        n_rows :: Int64
        hash_functions :: Vector{Function}
        tables :: Vector{DataStructures.DefaultDict{UInt, Vector{Int}}}
    end

    """
    Construct minhash index. 
    
    # Arguments
    - `threshold` the jaccard similarity above which documents are probably returned to be similar.
    - `shingle_size` the size of the shingles the strings will be cut into
    - `max_n_hashes` the maximum number of hash functions used
    """
    function MinHash(;threshold=0.5, shingle_size=3, max_n_hashes=200)
        n_bands, n_rows = _get_partition(threshold, max_n_hashes)
        n_hashes = n_bands * n_rows
        salts = [rand(UInt) for _ in 1:n_hashes]
        hash_functions = [str -> hash(str, salt) for salt in salts]
        tables = [DataStructures.DefaultDict{UInt, Vector{Int}}(Vector{Int}) for _ in 1:n_hashes]
        MinHash(threshold, Vector{String}(), shingle_size, n_bands, n_rows, hash_functions, tables)
    end

    function _get_partition(threshold, max_n_hashes)
        min_error = Inf
        opt = (0, 0)
        for b in 1 : max_n_hashes
            max_r = Int(floor(max_n_hashes / b))
            for r in 1 : max_r
            fp = _integrate( s -> 1 - (1 - s^r)^b, threshold, 0.0)
            fn = _integrate( s -> 1 - (1 - s^r)^b, threshold, 1.0)
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
    Push new string and save it's signature. Return calculated signature.
    """
    function push!(mh:: MinHash, s::AbstractString)::Vector{UInt}
        push!(mh.documents, s)
        signature = _calculate_signature(mh, s)
        index = length(mh.documents)
        for (signature_part, table) in zip(signature, mh.tables)
            # TODO: Optimizable, since this is already hashed?
            push!(table[signature_part], index)
        end
        signature
    end

    function _calculate_signature(mh:: MinHash, s::AbstractString)::Vector{UInt}
        shingles = _shingles(mh, s)
        [minimum((hash_func(s) for s in shingles)) for hash_func in mh.hash_functions]
    end

    function _shingles(mh, s)
        s_shingle = mh.shingle_size
        n = length(s)
        unique([s[thisind(s, i): thisind(s, j)] for (i,j) in 
            zip(1 : n - s_shingle + 1, s_shingle : n)])
    end

    """
    Find probably similar strings in index.
    """
    function find_similar(mh:: MinHash, s::AbstractString)::Vector{AbstractString}
        signature = _calculate_signature(mh, s)
        indices = vcat([table[signature_part] for (table, signature_part) in zip(mh.tables, signature)]...)
        mh.documents[unique(indices)]
    end

    """
    Find all pairs of similar strings in index.
    """
    function similar_pairs(mh:: MinHash)::Vector{Tuple{Int, Int}}
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

    export MinHash, push!, find_similar, similar_pairs
end # module
