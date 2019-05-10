using Revise
using LocalitySensitive
using Test

@testset "LocalitySensitive.jl" begin
    @testset "insert document find equal" begin
        lsh = LSH(shingle_size=3, n_hashes=200, threshold=0.9)
        insert!(lsh, "test")
        @test "test" == find_similar(lsh, "test")
        insert!(lsh, "58849")
        @test "58849" == find_similar(lsh, "58849")
    end

    @testset "insert document find any" begin
        lsh = LSH(shingle_size=3, n_hashes=200, threshold=0)
        insert!(lsh, "test")
        @test "test" == find_similar(lsh, "tests")
        insert!(lsh, "58849")
        @test "58849" == find_similar(lsh, "58849s")
    end
end
