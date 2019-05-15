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

