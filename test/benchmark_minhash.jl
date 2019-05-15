using LocalitySensitive
using BenchmarkTools

function fingerprint_all(docs::Vector{String}, mh::MinHash)
    for d in docs
       fingerprint(mh, shingle(d))
    end
end

@info "Benchmark shingling and fingerprinting"
lines = open("test/benchmark_minhash.jl") do f
    read(f, String)
end
documents = vcat([lines for _ in 1:200]...)
mh = MinHash()
@benchmark fingerprint_all(documents, mh)

@info "Benchmark push to index"
lines = open("test/benchmark_minhash.jl") do f
    read(f, String)
end
documents = vcat([lines for _ in 1:200]...)
mh = MinHash()
fingerprints = [fingerprint(mh, shingle(d)) for d in documents]
mhind = MinHashIndex(minhash=mh)
function pushall(fingerprints, mhind)
    for f in fingerprints
        push!(mhind, f)
    end
end
@benchmark pushall(fingerprints, mhind)

@info "Benchmark all similar pairs"
lines = open("test/benchmark_minhash.jl") do f
    readlines(f)
end
documents = vcat([lines for _ in 1:30]...)
mh = MinHash()
fingerprints = [fingerprint(mh, shingle(d)) for d in documents]
mhind = MinHashIndex(minhash=mh)
pushall(fingerprints, mhind)
@benchmark similar_pairs(mhind)

