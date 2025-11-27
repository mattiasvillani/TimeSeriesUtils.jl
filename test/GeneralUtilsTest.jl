# Testing funtions in GeneralUtils.jl
@testset "GeneralUtils.jl" begin
    @test nainterpret([1,missing,3])[2] == 2.0
    @test nainterpret!([1,missing,3])[2] == 2.0
end