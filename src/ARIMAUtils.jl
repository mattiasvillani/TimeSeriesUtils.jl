
""" 
    ARMAacf(; ϕ, θ, maxlag, pacf) 

Compute the Autocorrelation function (ACF) or partial ACF (PACF, if `pacf == true`) for an ARMA process for the lags k = 1,2,...,maxlag.

# Examples
```julia-repl
julia> acf = ARMAacf(ϕ = [0.8], θ = [0.0], maxlag = 10, pacf = false); 
julia> bar(0:10, acf)
```
""" 
function ARMAacf(; ϕ = [], θ = [], maxlag = maximum(length(ϕ), length(θ)+1), pacf = false)
    phi = ϕ; theta = θ;
    R"""
        ARMAacf = ARMAacf(ar = $phi, ma = $theta, lag.max = $maxlag, pacf = $pacf)
    """
    @rget ARMAacf
    if !pacf ARMAacf = ARMAacf[2:end] end
    return ARMAacf
end



""" 
    arma_reparam(x::Vector [, ztrans = "monahan", negative_signs = true]) 

Takes a p-dim vector `x` and returns parameters restricted to the have roots outside the unit circle, i.e. stationary(AR) or invertible (MA). The mapping is performed via the partial autocorrelations, P as: x -> P -> ϕ.

If `ztrans == "sigmoid"`, then the partial autocorrelations are parameterized as 
    P = (exp(x) - 1) / (exp(x) + 1 )

If `ztrans == "monahan"`, then the partial autocorrelations are parameterized as 
    P = x/√(1 + x²) 

If `negative_signs == true` 
coefficients are for polynomial 1 - ϕ₁B - ϕ₂B² - .... which is typically used for AR\\
If `negative_signs == false` 
coefficients are for polynomial
        1 + ϕ₁B + ϕ₂B² + .... which is typically used for MA

# Examples
```julia-repl
julia> ϕ, partials = arma_reparam([-4.0,4.0]; ztrans = "monahan");ϕ
2-element Vector{Float64}:
 -0.028966029557096595
 0.9701425001453319
julia> check_stationarity(ϕ)[1] # second element would return the eigenvalues
true
```
""" 
function arma_reparam(x; ztrans = "monahan", threshold = nothing, negative_signs = true)
    
    p = length(x)
    if ztrans == "sigmoid"
        P = (exp.(x) .- 1) ./ (exp.(x) .+ 1 )
        P[x .>= 700] = (1 .- exp.(-x[x .>= 700])) ./ (1 .+ exp.(-x[x .>= 700]))
    elseif ztrans == "monahan"
        P = x./sqrt.(1 .+ x.^2)
    elseif ztrans == "linear" # No transformation
        return x, NaN
    else
        error("ztrans must be either 'sigmoid' or 'monahan'")
    end

    if !isnothing(threshold)
        P = clamp.(P, -threshold, threshold)
    end
    
    if negative_signs 
        P = -P
    end
    ϕ = zeros(eltype(x), p, p) # Not sure we even need to allocate, but let's not worry about that.
    ϕ[1,1] = P[1]
    for k = 2:p
        for j = 1:k
            if k == j
                ϕ[k, j] = P[k]
            else
                ϕ[k, j] = ϕ[k-1, j] + P[k]*ϕ[k-1, k-j]
            end
        end
    end
    if negative_signs
        return -ϕ[end,:], P # returns ϕ in polynomial ϕ(B) = 1 - ϕ₁B - ϕ₂B² - .... used for AR
    else
        return ϕ[end,:], P # returns ϕ in polynomial ϕ(B) = 1 + ϕ₁B + ϕ₂B² + ... used for MA
    end
end

""" 
    inv_arma_reparam(ϕ; ztrans = "monahan") 

Converts from some AR parameters down to the unrestricted parameters. This inverts numerically, so inefficient.

""" 
function inv_arma_reparam(ϕ; ztrans = "monahan")

    f(θ,p) = arma_reparam(θ; ztrans)[1] - ϕ
    u0 = zeros(length(ϕ))
    prob = NonlinearProblem(f, u0)
    return solve(prob).u

end



""" 
    check_stationarity(ϕ::Vector) 

Test if the AR(p) polynomial 1 - ϕ₁B - ϕ₂B² - ... - ϕₚB^p corresponding to `ϕ` has all roots outside the unit circle, i.e. if the AR(p) process is stationary.

Returns a tuple where first element is `true` if stationary. Second element is a vector with the eigenvalues to the companion matrix.

# Examples
```julia-repl
julia> ϕ, P = arma_reparam([-4,4]; ztrans = "monahan") # returns ϕ which is stationary.
2-element Vector{Float64}:
 -0.028966029557096595
 0.9701425001453319
julia> check_stationarity(ϕ)[1]
true
```
""" 
function check_stationarity(ϕ::Vector)
    p = length(ϕ)
    companion = [ϕ' ; collect(I(p-1)) zeros(p-1,1)]
    eigen_companion = eigvals(companion)
    return all(abs.(eigen_companion) .< 1), eigen_companion
end

""" 
    sarma_reparam(θ::Vector, Θ::Vector, s, activeLags; ztrans = "monahan", threshold = nothing, negative_signs = true) 

Takes a p-dim vector `θ` with regular AR/MA coefficients and a P-dim vector with seasonal AR/MA coefficients `Θ` (with season `s`) and returns the *non-zero* coefficients in the product polynomial 
(1 - ϕ₁B - ϕ₂B² - ....)(1 - Φ₁B - Φ₂B² - ....) = (1 - ψ₁B - Ψ₂B² - ....) 
where both sets of parameters (ϕ and Φ) are restricted to the have roots outside the unit circle, i.e. stationary(AR) or invertible (MA). The mapping is performed via the partial autocorrelations, P as: θ -> P₁ -> ϕ and Θ -> P₂ -> Φ. 

If `ztrans == "sigmoid"`, then the partial autocorrelations are parameterized as 
    P = (exp(x) - 1) / (exp(x) + 1 )

If `ztrans == "monahan"`, then the partial autocorrelations are parameterized as 
    P = x/√(1 + x²) 

If `negative_signs == true` 
coefficients are for polynomial 1 - ϕ₁B - ϕ₂B² - .... which is typically used for AR\\
If `negative_signs == false` 
coefficients are for polynomial
        1 + ϕ₁B + ϕ₂B² + .... which is typically used for MA

""" 
function sarma_reparam(θ, Θ, s, activeLags = nothing; ztrans = "monahan", 
        threshold = nothing, negative_signs = true)

    p = length(θ)
    P = length(Θ)
    if isnothing(activeLags)
        activeLags = FindActiveLagsSAR(p, P, s)
    end
    ϕ = arma_reparam(θ; ztrans = ztrans, threshold = threshold, 
        negative_signs = negative_signs)[1] .+ eps()
    Φ = arma_reparam(Θ; ztrans = ztrans, threshold = threshold, 
        negative_signs = negative_signs)[1] .+ eps()
    ARpoly =  Polynomial([1;-ϕ], :z)
    ARseasonpolymat = [zeros(s-1,P);-Φ']; 
    ARseasonpoly =  Polynomial([1;ARseasonpolymat[:]], :z)
    ϕ̃ = -coeffs(ARpoly*ARseasonpoly)[2:end]

    return ϕ̃[activeLags]

end


# Compute AR parameters from partial autocorrelations
function arma_reparam_partials(P::Vector; negative_signs = true)
    p = length(P)
    if negative_signs
        P = -P
    end
    ϕ = zeros(eltype(P), p, p) # Not sure we even need to allocate, but let's not worry about that.
    ϕ[1,1] = P[1]
    for k = 2:p
        for j = 1:k
            if k == j
                ϕ[k, j] = P[k]
            else
                ϕ[k, j] = ϕ[k-1, j] + P[k]*ϕ[k-1, k-j]
            end
        end
    end
    if negative_signs
        return -ϕ[end,:] # returns ϕ in polynomial ϕ(B) = 1 - ϕ₁B - ϕ₂B² - .... used for AR
    else
        return ϕ[end,:] # returns ϕ in polynomial ϕ(B) = 1 + ϕ₁B + ϕ₂B² + ... used for MA
    end
end


""" 
    Arima(y; order = [0,0,0], seasonal = [0,0,0], xreg = nothing, include_mean = true,
        include_drift = false, include_constant = true, frequency = 1, deltat = 1) 

R function that estimates an ARIMA model for the time series `y`.
""" 
function Arima(y; order = [0,0,0], seasonal = [0,0,0], xreg = nothing, include_mean = true,
    include_drift = false, include_constant = true, frequency = 1, deltat = 1)
    R"""
    suppressMessages(library(forecast))
    fittedModel = Arima(ts($y, frequency = $frequency, deltat = $deltat), order = $order, seasonal = $seasonal, xreg = $xreg, 
        include.mean = $include_mean, include.drift = $include_drift, include.constant = $include_constant)
    """
    @rget fittedModel
    return fittedModel
end

""" 
    simARMA(ϕ, θ, c, σ, T)

Simulate a time series of length `T` from the univariate ARMA process with the vector of AR coefficients `ϕ`, vector of MA coefficients `θ`, constant (intercept) `c` and innovation standard deviation `σ`: 

```math
x_t = c + \\sum_{j=1}^p \\phi_j x_{t-j} + \\sum_{k=1}^q \\theta_k a_{t-k} 
``` 

Uses the [ARCHModels.jl](https://github.com/s-broda/ARCHModels.jl) package.

# Examples
Simulate `T=5` observations from the ARMA(1,1) process:
```julia-repl
julia> simARMA([0.7], [0.3], 0.0, 1.0, 5)
5-element Vector{Float64}:
  0.30639811188691285
  0.8315918084473066
  1.1336754037240127
 -0.15634038037829656
 -0.3521705897607527
```
""" 
function simARMA(ϕ, θ, c, σ, T)
    dgp = UnivariateARCHModel(
        ARCH{0}([σ^2]), 
        randn(T); 
        dist = StdNormal(), 
        meanspec = ARMA{length(ϕ),length(θ)}(vcat(c, ϕ, θ))
    )
    simulate!(dgp; warmup = Int(ceil(0.1*T)))
    return dgp.data
end

""" 
    ℓARMA(data::Vector, ϕ, θ, c, σ)

Compute the time domain likelihood of the `data` in the univariate ARMA process with the vector of AR coefficients `ϕ`, vector of MA coefficients `θ`, constant (intercept) `c` and innovation standard deviation `σ`: 

```math
x_t = c + \\sum_{j=1}^p \\phi_j x_{t-j} + \\sum_{k=1}^q \\theta_k a_{t-k} 
```

Uses the [ARCHModels.jl](https://github.com/s-broda/ARCHModels.jl) package.

# Examples
log-likelihood for a time series with 5 observeration in the AR(1) process with ϕ = 0.8, c = 0, σ = 1:
```julia-repl
julia> ℓARMA([1.1,2.2,1.7,1.9,2.3], [0.8], [0.0], 0.0, 1.0)
-6.191492666023364
```
"""
function ℓARMA(data::Vector, ϕ, θ, c, σ)
    p = length(ϕ)
    q = length(θ)
    armaModel = UnivariateARCHModel(
        ARCH{0}([σ^2]), 
        data; 
        dist = StdNormal(), 
        meanspec = ARMA{p,q}(vcat(c, ϕ, θ))
    )
    return loglikelihood(armaModel)
end

""" 
    SpecDensARMA(ω, ϕ, θ, σ²) 

Compute spectral density for the univariate ARMA model over domain ω ∈ [-π,π]. 

- ω is a radial frequency
- ϕ is a vector of AR coefficients
- θ is a vector of MA coefficients
- σ² is the noise variance

# Examples
The spectral density for an AR(1) process with unit noise variance is
```doctests 
julia> SpecDensARMA(0.5, 0.9, 0, 1)
0.6909224383713601
```
""" 
function SpecDensARMA(ω, ϕ, θ, σ²)
    ARpoly =  Polynomial([1;-ϕ], :z)
    MApoly =  Polynomial([1;θ], :z) 
    specDens = (σ²/(2π))*(abs(MApoly(exp(-im*ω)))^2/abs(ARpoly(exp(-im*ω)))^2)
	return specDens
end 

function SpecDensARMA(ω::AbstractVector, ϕ, θ, σ²)
    ARpoly =  Polynomial([1;-ϕ], :z)
    MApoly =  Polynomial([1;θ], :z) 
    specDens = (σ²/(2π))*(abs.(MApoly.(exp.(-im*ω))).^2 ./ abs.(ARpoly.(exp.(-im*ω))).^2)
	return specDens
end 

""" 
    SpecDensSARMA(ω, ϕ, θ, Φ, Θ, σ², s) 

Compute spectral density for the seasonal univariate ARMA model over domain ω ∈ [-π,π]. 

- ω is a radial frequency
- ϕ is a vector of AR coefficients
- θ is a vector of MA coefficients
- Φ is a vector of seasonal AR coefficients
- Θ is a vector of seasonal MA coefficients
- σ² is the noise variance
- s is the seasonality (e.g. s = 12 for monthly data)

# Examples
The spectral density for an AR(2)(1)₄ process with unit noise variance is
```doctests 
julia> SpecDensSARMA(0.5, [0.8,-0.2], 0, 0.6, 0, 1, 4)
0.4053556832813598
```
""" 
function SpecDensSARMA(ω, ϕ, θ, Φ, Θ, σ², s)
    ARpoly =  Polynomial([1;-ϕ], :z)
    MApoly =  Polynomial([1;θ], :z) 
    P = length(Φ)
    Q = length(Θ)
    ARpolymat = [zeros(s-1,P);-Φ']; 
    MApolymat = [zeros(s-1,Q);Θ'];
    ARseasonpoly =  Polynomial([1;ARpolymat[:]], :z)
    MAseasonpoly =  Polynomial([1;MApolymat[:]], :z)
    specDens = (σ²/(2π))*(abs(MApoly(exp(-im*ω)))^2/abs(ARpoly(exp(-im*ω)))^2)* 
        (abs(MAseasonpoly(exp(-im*ω)))^2/abs(ARseasonpoly(exp(-im*ω)))^2)
	return specDens
end 

function SpecDensSARMA(ω::AbstractVector, ϕ, θ, Φ, Θ, σ², s)
    ARpoly =  Polynomial([1;-ϕ], :z)
    MApoly =  Polynomial([1;θ], :z) 
    P = length(Φ)
    Q = length(Θ)
    ARpolymat = [zeros(s-1,P);-Φ']; 
    MApolymat = [zeros(s-1,Q);Θ'];
    ARseasonpoly =  Polynomial([1;ARpolymat[:]], :z)
    MAseasonpoly =  Polynomial([1;MApolymat[:]], :z)
    specDens =(σ²/(2π))*(abs.(MApoly.(exp.(-im*ω))).^2 ./ abs.(ARpoly.(exp.(-im*ω))).^2) .* 
        (abs.(MAseasonpoly.(exp.(-im*ω))).^2 ./ abs.(ARseasonpoly.(exp.(-im*ω))).^2)
	return specDens
end 

""" 
    SpecDensMultiSARMA(ω, ϕ, θ, σ², s) 

Compute spectral density for the multi-seasonal univariate ARMA model over domain ω ∈ [-π,π]. ω can be a scalar or a vector.

- ω is a radial frequency
- ϕ is a vector of of vectors where ϕ[l] are the AR coefficients in l:th polynomial
- θ is a vector of of vectors where θ[l] are the MA coefficients in l:th polynomial
   empty vector in ϕ means that that specific polynomial is not present in AR. Same for θ.
- σ² is the noise variance
- s is the seasonality (e.g. s = (1,4,12) for regular, quarterly and monthly seasons)

# Examples
The spectral density at ω = 0.5 for an SAR(p = (2,1), q = (0,0),s=(1,4)) process with unit noise variance is
```doctests 
julia> SpecDensMultiSARMA(0.5, [[0.8,-0.2],[0.6]], [[],[]], 1, [2,1], [0,0], [1 4])
```
""" 
function SpecDensMultiSARMA(ω, ϕ, θ, σ², s)
    
    
    nSeasons = length(s)
    p = length.(ϕ)
    q = length.(θ)
    length(p) == length(q) == length(s) ? nothing : error("length(ϕ) and length(θ) must equal lenght(s). Use empty elements if AR or MA lacks one or more seasonal polynomials")

    ARpolynomials = Array{Polynomial}(undef, nSeasons)
    count = 0
    for l in 1:nSeasons
        ARseasonpolymat = [zeros(s[l]-1, p[l]);-ϕ[l]']; 
        ARpolynomials[l] = Polynomial([1;ARseasonpolymat[:]], :z)
        count = count + p[l]
    end

    
    MApolynomials = Array{Polynomial}(undef, nSeasons)
    count = 0
    for l in 1:nSeasons
        MAseasonpolymat = [zeros(s[l]-1, q[l]);θ[l]']; 
        MApolynomials[l] = Polynomial([1;MAseasonpolymat[:]], :z)
        count = count + q[l]
    end

    specDens = (σ²/(2π))*abs(prod(MApolynomials)(exp(-im*ω)))^2 / 
        abs(prod(ARpolynomials)(exp(-im*ω)))^2
	return specDens

end 

function SpecDensMultiSARMA(ω::AbstractVector, ϕ, θ, σ², s)

    nSeasons = length(s)
    p = length.(ϕ)
    q = length.(θ)
    length(p) == length(q) == length(s) ? nothing : error("length(ϕ) and length(θ) must equal lenght(s). Use empty elements if AR or MA lacks one or more seasonal polynomials")

    ARpolynomials = Array{Polynomial}(undef, nSeasons)
    count = 0
    for l in 1:nSeasons
        ARseasonpolymat = [zeros(s[l]-1, p[l]);-ϕ[l]']; 
        ARpolynomials[l] = Polynomial([1;ARseasonpolymat[:]], :z)
        count = count + p[l]
    end

    MApolynomials = Array{Polynomial}(undef, nSeasons)
    count = 0
    for l in 1:nSeasons
        MAseasonpolymat = [zeros(s[l]-1, q[l]);θ[l]']; 
        MApolynomials[l] = Polynomial([1;MAseasonpolymat[:]], :z)
        count = count + q[l]
    end

    specDens = (σ²/(2π))*abs.(prod(MApolynomials).(exp.(-im*ω))).^2 ./ 
        abs.(prod(ARpolynomials).(exp.(-im*ω))).^2

	return specDens
end 


""" 
    SpecDensARTFIMA(ω, ϕ, θ, d, λ, σ²) 

Compute spectral density for the univariate ARTFIMA model over domain ω ∈ [-π,π]. 

- ω is a radial frequency
- ϕ is a vector of AR coefficients
- θ is a vector of MA coefficients
- d is the fractional differenting parameter
- λ ≥ 0 is the tempering parameter 
- σ² is the noise variance

# Examples
The spectral density for an AR(1) process with unit noise variance is
```doctests 
julia> SpecDensARTFIMA(0.5, 0.9, 0, 0, 0, 1)
0.6909224383713601
```
""" 
function SpecDensARTFIMA(ω, ϕ, θ, d, λ, σ²)
    ARpoly =  Polynomial([1;-ϕ], :z)
    MApoly =  Polynomial([1;θ], :z) 
    specDens = (σ²/(2π))*(abs(MApoly(exp(-im*ω)))^2/abs(ARpoly(exp(-im*ω)))^2)*
		abs(1-exp(-(λ+im*ω)))^(-2*d)
	return specDens
end 


"""
    Simulate from a uniform prior for AR(p) over stationary region by simulating partials.
"""
function sim_uniformAR(p, nsim, trans = "monahan")
    if (trans != "monahan")
        error("Only 'monahan' transformation is currently implemented")
    end     

    Psim = zeros(nsim,p)
    ϕsim = zeros(nsim,p)
    betadists = []
    for k = 1:p
        push!(betadists, -1 + 2*Beta((k+1)/2, floor(k/2) + 1))
    end
    for i = 1:nsim
        Psim[i,:] = rand.(betadists)
        ϕsim[i,:] = arma_reparam_partials(Psim[i,:])
    end
    return ϕsim, Psim
end