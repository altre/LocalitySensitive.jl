abstract type LSHIndex end

function _get_partition(threshold, max_n_hashes)
    min_error = Inf
    opt = (0, 0)
    bestfp, bestfn = (0,0)
    for b in 1 : max_n_hashes
        max_r = Int(floor(max_n_hashes / b))
        for r in 1 : max_r
            fp = _integrate(s -> 1 - (1 - s^r)^b, 0.0, threshold)
            fn = _integrate(s -> (1 - s^r)^b, threshold, 1.0)
            error = fp + fn
            if error < min_error
                min_error = error
                opt = (b, r)
                bestfp, bestfn = fp, fn
            end
        end
    end
    @show bestfp, bestfn
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
Find all pairs of similar strings in index.
"""
function similar_pairs(mh::LSHIndex)::Vector{Tuple{Int, Int}}
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
