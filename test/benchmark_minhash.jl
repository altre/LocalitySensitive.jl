using LocalitySensitive
using BenchmarkTools

benchmark_documents = readlines("../resources/benchmark_data.csv")

function fingerprint_benchmark(docs::Vector{String}, mh::MinHash)
    for d in docs
       fingerprint(mh, shingle(d))
    end
end

@info "Benchmark shingling and fingerprinting"
mh = MinHash()
@benchmark fingerprint_benchmark(benchmark_documents, mh)
shingles = [shingle(d) for d in benchmark_documents]
@benchmark fingerprint_all(mh, shingles)

@info "Benchmark push to index"
fingerprints = [fingerprint(mh, shingle(d)) for d in benchmark_documents]
mhind = MinHashIndex(minhash=mh)
function pushall(fingerprints, mhind)
    for f in fingerprints
        push!(mhind, f)
    end
end
@benchmark pushall(fingerprints, mhind)

@info "Benchmark all similar pairs"
mhind = MinHashIndex(minhash=mh)
pushall(fingerprints, mhind)
@benchmark similar_pairs(mhind)

pairs = similar_pairs(mhind)