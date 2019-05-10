using Revise
using LocalitySensitive
using Test

@testset "LocalitySensitive.jl" begin
    @testset "push document find equal" begin
        mh = MinHash(threshold=0.9)
        push!(mh, "test")
        @test ["test"] == find_similar(mh, "test")
        push!(mh, "58849")
        @test ["58849"] == find_similar(mh, "58849")
    end

    @testset "push document find any" begin
        mh = MinHash(threshold=0)
        push!(mh, "test")
        @test ["test"] == find_similar(mh, "testsim")
        push!(mh, "58849")
        @test ["58849"] == find_similar(mh, "58849sim")
    end
    
    @testset "Find all pair indices" begin
        mh = MinHash(threshold=0)
        push!(mh, "test")
        push!(mh, "testSIM")
        push!(mh, "SIMtest")
        push!(mh, "58849")
        push!(mh, "58849sim")
        push!(mh, "sim58849")
        @test sort([(1,2),(1,3),(2,3),(4,5),(4,6),(5,6)]) == sort(similar_pairs(mh))
    end
end
