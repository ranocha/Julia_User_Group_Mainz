using Test
using MyPkg

@testset "Testing foo2" begin
    @test foo2(2) == 8
    @test foo2(3) == 27
    @test foo2(4) == 64
    @test foo2(5) == 20
end
