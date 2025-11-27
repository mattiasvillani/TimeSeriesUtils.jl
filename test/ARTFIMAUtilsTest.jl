@testset "ARTFIMAUtils.jl" begin
    # White noise should have constant density of 1/(2π)
    σ² = 2.1
    ω₁ = π/2
    ω₂ = π/4
    s = 12
    @test SpecDensARTFIMA(ω₁, 0.0, 0, 0, 0, σ²) == SpecDensARTFIMA(ω₂, 0.0, 0, 0, 0, σ²) ==
        σ²/(2π)

    # ARIMA is a special case of ARTFIMA
    @test SpecDensARMA(0.42, [0.5,-0.2], [0.5,-0.2], 0.25) == 
        SpecDensARTFIMA(0.42, [0.5,-0.2], [0.5,-0.2], 0, 0, 0.25)
    
end