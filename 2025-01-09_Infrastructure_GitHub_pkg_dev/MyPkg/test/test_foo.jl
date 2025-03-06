using Test
using MyPkg

@testset "Testing foo" begin
    @test foo(2) == 4
    @test foo(3) == 9
    @test foo(4) == 16
    @test foo(5) == 25
end
