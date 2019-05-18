using LocalitySensitive
using Test
using SparseArrays

@testset "SimHash" begin
    @testset "Estimate cosine distance" begin
        sh = SimHash(1)
        f = fingerprint(sh, [("ab", 1), ("bc",1)])
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
        docs = ["test","testSIM","SIMtest","58849","58849sim","sim58849"]
        docs = [[(s, 1) for s in shingle(d)] for d in docs]
        for doc in docs
            push!(shind, fingerprint(sh, doc))
        end
        @test sort([(1,2),(1,3),(2,3),(4,5),(4,6),(5,6)]) == sort(similar_pairs(shind))
    end
end