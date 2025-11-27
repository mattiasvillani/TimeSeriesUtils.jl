module TimeSeriesUtils

using Plots, Distributions, FFTW, Polynomials
using LaTeXStrings
using RCall

include("GeneralUtils.jl")
export ts, stl, mstl, nainterpret!, nainterpret

include("SpectralUtils.jl")
export spectrum, pspectrum, pgram, pGramDistribution
export ℓwhittle
export simProcessSpectral

include("ARIMAUtils.jl")
export ARMAacf, arma_reparam, inv_arma_reparam, check_stationarity
export sarma_reparam
export Arima, simARMA, ℓARMA
export SpecDensARMA, SpecDensSARMA, SpecDensMultiSARMA

include("ARTFIMAUtils.jl")
export artfima, artsim, artfima_pred, SpecDensARTFIMA

end
