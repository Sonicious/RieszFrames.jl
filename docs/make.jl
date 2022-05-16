# make src directory accessible
#push!(LOAD_PATH, "../src/")

using Documenter
using RieszFrames

makedocs(
    sitename="RieszFrames",
    format=Documenter.HTML(),
    modules=[RieszFrames],
    pages=[
        "Introduction" => "index.md",
        "Algorithms" => "algorithms.md",
        "Utility Functions" => "utilities.md",
        "Frame Functions" => "frames.md"
    ]
)

deploydocs(
    devbranch="main",
    repo="git@github.com:Sonicious/RieszFrames.jl.git"
)
