@testset "ARIMAUtils.jl" begin
    # White noise should have constant density of σ²/(2π)
    σ² = 2.1
    ω₁ = π/2
    ω₂ = π/4
    s = 12
    @test SpecDensARMA(ω₁, 0, 0, σ²) == SpecDensARMA(ω₂, 0, 0, σ²) == σ²/(2π)
    @test SpecDensSARMA(ω₁, 0, 0, 0, 0, σ², s) == SpecDensSARMA(ω₂, 0, 0, 0, 0, σ², s) == 
        σ²/(2π)
    @test SpecDensMultiSARMA(ω₁, 0, 0, σ², s) == SpecDensMultiSARMA(ω₂, 0, 0, σ², s) == 
        σ²/(2π)

end