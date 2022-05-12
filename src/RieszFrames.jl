module RieszFrames

# deps
using FFTW

# global 
global const rieszbank = Ref{Vector{Array}}()
global const framebank = Ref{Vector{Array}}()

#files and exports
include("FrameFunctions.jl")
export Held, Papadakis, Shannon, Simoncelli

include("FilterGeneration.jl")
export GenerateRadius2D, GenerateRadius1D, GenerateRieszFilters2D, GenerateFrameFilters2D

include("EoB.jl")
export Eob, EobAlgorithm, EobInitialization

end #module
