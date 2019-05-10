module LocalitySensitive
    import Parameters
    import Base.insert!

    Parameters.@with_kw struct LSH
        shingle_size::Int = 3
        n_hashes::Int = 200
        threshold::Float64 = 0.5
    end
    
    function insert!(lsh:: LSH, s::AbstractString)

    end

    function find_similar(lsh:: LSH, s::AbstractString)
        return s
    end

    export LSH, insert!, find_similar
end # module
