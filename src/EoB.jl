"""
    Eob(input[, type=Float64])

Equalization of Brighness (Eob) with default values and new Initializationof the steerable filterbanks.
# Arguments
- `input`: this is the input grayscale image of any size. It will be converted into an according floating point structure during runtime
- `type`: The type specifies the floating point type of the filterbanks (default: Float64).
# Example
```jldoctest
julia> Eob(rand(Float32,16,16));
```
"""
function Eob(input, type=Float64)
    maxscale = 8
    framefunction = Held
    regularization = 0.001
    try
        input = convert(Array{type,2}, input)
    catch
        throw(DomainError(input, "Please use data which can be converted to a 2D float Array"))
    end
    EobInitialization(size(input), maxscale, type, framefunction)
    return EobAlgorithm(input, regularization, framebank[], rieszbank[], type)
end

"""
    EobInitialization(inputsize[, maxscale=8, type=Float64, framefunction=Held])

Calculate the Filters and all preallocations for the EoB Algorithm for images of size in input

This function is especially usefulf for the generation of a new filterbank which is then used for batch processing of the Eob Algorithm. All important parameters can be set individually. The function wraps the calls to `GenerateFrameFilters2D`
The functions of [Generation Functions](@ref) give more explanation about the parameters
# Example
```jldoctest EobTest
julia> testImage = rand(Float64,16,16);
julia> EobInitialize(size(testImage));
```
"""
function EobInitialization(inputsize, maxscale=8, type=Float64, framefunction=Held)
    framebank[] = GenerateFrameFilters2D(inputsize, maxscale, framefunction)
    rieszbank[] = GenerateRieszFilters2D(inputsize, 1)
    nothing
end

"""
    EobAlgorithm(input[, regularization=0.001, framefilterbank=framebank[], rieszfilterbank=rieszbank[], type=Float64])

Calculate the EqualizationOfBrightness Algorithm with given values. It uses the default Values from the Initialization

This function is for more sophisticated applications, where the algorithm is batch-processed and the filterbanks are reused every time. This improves the performance because they can be reused. additionally the regularization parameter can be set individually. The applicaiton of the filterbanks is multithreaded
# Arguments
- `input`: input image. Must be a 2D Array of a floating point type
- `regularization`: this parameter controls the regularization of the algorithm. 0 is a possible input for no regularization (default: 0.001)
- `framefilterbank`: a valid frame filterbank array for the input image (default: global initialized framefilterbank)
- `rieszfilterbank`: a valid riesz filterbank array for the input image (default: global initialized rieszfilterbank)
- `type`: The type specifies the floating point type of the result (default: Float64).
# Example
```jldoctest EobTest
julia> result = EobAlgorithm(testImage);
```
"""
function EobAlgorithm(input, regularization=0.001, framefilterbank=framebank[], rieszfilterbank=rieszbank[], type=Float64)
    result = zeros(complex(type), size(input))
    # calc the fft first and shift to fit the filters
    spectrum = fftshift(fft(input))
    # apply the filterbank to the spectrum
    Threads.@threads for filter_idx = 1:length(framefilterbank)
        nonriesz = spectrum .* framefilterbank[filter_idx]
        riesz1 = real(ifft(fftshift(nonriesz .* rieszfilterbank[1])))
        riesz2 = real(ifft(fftshift(nonriesz .* rieszfilterbank[2])))
        nonriesz = real(ifft(fftshift(nonriesz)))
        normal = hypot.(riesz1, riesz2)
        amplitude = hypot.(nonriesz, normal)
        amplitude .= amplitude .* amplitude
        cosine = atan.(normal, nonriesz)
        regularized = cos.(cosine) .* (amplitude) ./ (amplitude .+ regularization)
        result .+= framefilterbank[filter_idx] .* fftshift(fft(regularized))
    end
    result = real(ifft(fftshift(result)))
    return result
end