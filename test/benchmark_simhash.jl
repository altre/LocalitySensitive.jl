using LocalitySensitive
using BenchmarkTools
using TextAnalysis
using SparseArrays

benchmark_documents = readlines("resources/benchmark_data.csv")

@info "Benchmark fingerprinting"
function fingerprint_benchmark(docs, sh::SimHash)
    for d in docs
       fingerprint(sh, d)
    end
end
sh = SimHash()
els = [[(s,1) for s in shingle(d)] for d in benchmark_documents]
@benchmark fingerprint_benchmark(els, sh)

names_doc = [StringDocument(doc) for doc in benchmark_documents]
crps = Corpus(names_doc)
update_lexicon!(crps)
document_matrix = DocumentTermMatrix(crps)
tfidf = SparseMatrixCSC(tf_idf(document_matrix)')

@benchmark fingerprint_all(sh, tfidf)

fps = fingerprint_all(sh, tfidf)


@info "Push to index"
function pushall(fingerp)
    sh = SimHashIndex(0.85)
    for f in fingerp
        push!(sh, f)
    end
    sh
end

# @benchmark pushall(fp)
shind = pushall(fps)
estimate_cosine(fps[24], fps[25])
find_similar(shind, fps[24])
benchmark_documents[find_similar(shind, fps[24])]
benchmark_documents[find_similar(shind, fps[1])]
benchmark_documents[1]


function simall(fingerp, sh::SimHashIndex)
    for f in fingerp
        find_similar(sh, f)
    end
end
using Profile
Profile.clear()
Profile.init(n=10000000, delay=0.00001)
@benchmark simall(fp, sh)
@profile simall(fp[1:3000], sh)
Profile.print()
using ProfileView
ProfileView.view()
pairs = similar_pairs(sh)