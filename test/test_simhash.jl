using LocalitySensitive
using Test
using SparseArrays
using TextAnalysis

@testset "SimHash" begin
    @testset "Estimate cosine similarity" begin
        sh = SimHash(1)
        features_weights =  [("ab", 1), ("bc",1)]
        f = fingerprint(sh, features_weights)
        @test 1 == estimate_cosine(f, f)
        f2 = fingerprint(sh, [("fg", 1), ("de",1)])
        @test estimate_cosine(f, f2) < 0.2
    end
    
    @testset "Fingerprint sparse vector" begin
        sh = SimHash(1)
        f = fingerprint(sh, sparse([1, 0]))
        @test 1 == estimate_cosine(f, f)
        f2 = fingerprint(sh, sparse([0, 1]))
        @test estimate_cosine(f, f2) < 0.1
    end

    @testset "Find similar" begin
        sh = SimHash(1)
        shind = SimHashIndex(0.2)
        f = fingerprint(sh, sparse([1, 0]))
        push!(shind, f)
        @test [] == find_similar(shind, fingerprint(sh, sparse([0, 1])))
        @test [1] == find_similar(shind, fingerprint(sh, sparse([1, 0])))
        @test [1] == find_similar(shind, fingerprint(sh, sparse([0.5, 0])))
    end
    
    @testset  "Find all pair indices" begin
        sh = SimHash(1)
        shind = SimHashIndex(0.3)
        @test [] == similar_pairs(shind)
        docs = [
            "this is an english sentence", 
            "this is an english sentence this is an english sentence", 
            "this too, is an english sentence", 
            "no resemblance whatsoever"
        ]
        docs = Document.(docs)
        crps = Corpus(docs)
        update_lexicon!(crps)
        m = DocumentTermMatrix(crps)
        features = tf_idf(m)
        for i in 1:length(docs)
            push!(shind, fingerprint(sh, features[i,:]))
        end
        @test [(1, 2), (1, 3), (2, 3)] == sort(similar_pairs(shind))
    end
end