global rieszbank = Ref{Vector{Array}}
global framebank = Ref{Vector{Array}}

"""
Eob

Equalization of Brighness (Eob) with default values
"""
function Eob(input)

end

"""
EobInitialization

Calculate the Filters and all preallocations for the EoB Algorithm for images of size in input
"""
function EobInitialization(input)
    inputsize = size(input)
end

"""
EobAlgorithm

Calculate the EqualizationOfBrightness Algorithm with given values. It uses the default Values from hte Initialization
"""
function EobAlgorithm(input, framefilterbank=rieszbank, framefilterbank=framebank)

end