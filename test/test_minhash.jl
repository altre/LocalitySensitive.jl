using LocalitySensitive
using Test
using Random

@testset "MinHash" begin
    @testset "estimate jaccard distance" begin
        mh = MinHash()
        test_fp = fingerprint(mh, shingle("test", size=3))
        @test estimate_jaccard(test_fp, test_fp) == 1
        @test estimate_jaccard(test_fp, fingerprint(mh, shingle("1234", size=3))) == 0
        @test estimate_jaccard(test_fp, fingerprint(mh, shingle("testr", size=3))) > 0.5

        @test estimate_jaccard(fingerprint(mh, [1,2,3,4]), fingerprint(mh, [1,2,3,4])) == 1
        @test estimate_jaccard(fingerprint(mh, [1,2,3,4]), fingerprint(mh, [1,2,3])) > 0.5
        @test estimate_jaccard(fingerprint(mh, [1,2,3,4]), fingerprint(mh, [5,6])) == 0
        
        @test estimate_jaccard(fingerprint(mh, 1:4), fingerprint(mh, 1:4)) == 1
        @test estimate_jaccard(fingerprint(mh, 1:4), fingerprint(mh, 5:6)) == 0

        sets = [collect(i:j) for (i,j) in zip(1:100, 3:102)]
        rng = MersenneTwister(1)
        for i in 1:10
            s1 = sets[randperm(rng, length(sets))][1:i*10]
            @test (i-1)/10 < estimate_jaccard(fingerprint(mh, sets), fingerprint(mh, s1)) < (i+1)/10
        end
    end

    @testset "Special function for strings" begin
        mh = MinHash()
        test_fp = fingerprint(mh, "test")
        @test estimate_jaccard(test_fp, test_fp) == 1
        @test estimate_jaccard(test_fp, fingerprint(mh, "1234")) == 0
        @test estimate_jaccard(test_fp, fingerprint(mh, "testr")) > 0.5
    end

    @testset "push document to index" begin
        mh = MinHash()
        mhind = MinHashIndex(fingerprinter=mh, threshold=0.1)
        signature = fingerprint(mh, "test")
        push!(mhind, signature)
        signature = fingerprint(mh, "testsim")
        @test [1] == find_similar(mhind, signature)
        push!(mhind, fingerprint(mh,"58849"))
        signature = fingerprint(mh, "58849sim")
        @test [2] == find_similar(mhind, signature)
        
        sets = [collect(i:j) for (i,j) in zip(1:100, 3:102)]
        mh = MinHash(400, MersenneTwister(1))
        for i in 3:6
            mhind = MinHashIndex(fingerprinter=mh, threshold=(i-2)/10)
            push!(mhind, fingerprint(mh, sets))
            s1 = sets[1:i*10]
            @test [1] == find_similar(mhind, fingerprint(mh, s1))
            mhind = MinHashIndex(fingerprinter=mh, threshold=(i+3)/10)
            push!(mhind, fingerprint(mh, sets))
            @test [] == find_similar(mhind, fingerprint(mh, s1))
        end
    end
    
    @testset "Find all pair indices" begin
        mh = MinHash()
        mhind = MinHashIndex(fingerprinter= mh, threshold=0)
        push!(mhind, fingerprint(mh, "test"))
        push!(mhind, fingerprint(mh, "testSIM"))
        push!(mhind, fingerprint(mh, "SIMtest"))
        push!(mhind, fingerprint(mh, "58849"))
        push!(mhind, fingerprint(mh, "58849sim"))
        push!(mhind, fingerprint(mh, "sim58849"))
        @test sort([(1,2),(1,3),(2,3),(4,5),(4,6),(5,6)]) == sort(similar_pairs(mhind))
    end

    @testset "Fingerprint all" begin
        mh = MinHash()
        documents = ["als", "kjsqer", "al;kjslkj", "sld"]
        fingerprints = fingerprint_all(mh, [shingle(d, size=3) for d in documents])
        mhind = MinHashIndex(fingerprinter=mh, threshold = 0.9)
        for f in fingerprints
            push!(mhind, f)
        end
        @test [1] == find_similar(mhind, fingerprint(mh, "als"))
        @test [] == find_similar(mhind, fingerprint(mh, "al"))
    end
end
