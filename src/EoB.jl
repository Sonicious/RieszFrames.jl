"""
Eob

Equalization of Brighness (Eob) with default values
"""
function Eob(input, type=Float64)
    try
        input = convert(Array{type,2}, input)
    catch
        throw(DomainError(input, "Please use data which can be converted to a 2D float Array"))
    end
    EobInitialization(size(input))
    return EobAlgorithm(input)
end

"""
EobInitialization

Calculate the Filters and all preallocations for the EoB Algorithm for images of size in input
"""
function EobInitialization(inputsize, maxscale=8, type=Float64, framefunction=Held)
    framebank[] = GenerateFrameFilters2D(inputsize, maxscale, framefunction)
    rieszbank[] = GenerateRieszFilters2D(inputsize, 1)
    fftforward[] = plan_fft(rand(complex(type), inputsize), 2, flags=FFTW.MEASURE)
    fftinverse[] = plan_ifft(rand(complex(type), inputsize), 2, flags=FFTW.MEASURE)
    nothing
end

"""
EobAlgorithm

Calculate the EqualizationOfBrightness Algorithm with given values. It uses the default Values from hte Initialization
"""
function EobAlgorithm(input, regularization = 0.001, framefilterbank=framebank[], rieszfilterbank=rieszbank[], type=Float64)
    result = zeros(size(input))
    # calc the fft first and shift to fit the filters
    spectrum = fftshift(fft(input))
    # apply the filterbank to the spectrum
    for filter_idx = 1:length(framefilterbank)
        nonriesz = spectrum .* framefilterbank[filter_idx]
        riesz1 = real(ifft(fftshift(nonriesz .* rieszfilterbank[1])))
        riesz2 = real(ifft(fftshift(nonriesz .* rieszfilterbank[2])))
        nonriesz = real(ifft(fftshift(nonriesz)))
        normal = hypot.(riesz1, riesz2)
        amplitude = hypot.(nonriesz, normal)
        amplitude = amplitude .* amplitude
        cosine = atan.(normal, nonriesz)
        regularized = cos.(cosine) .* (amplitude) ./ (amplitude .+ regularization)
        result += framefilterbank[filter_idx] .* fftshift(fft(regularized))
    end
    result = real(ifft(fftshift(result)))
    return result
end