"""
GenerateRieszFilters

Generate the Frame filterbank for the according filtersize and choosen Motherwavelet
"""
function generateFrameFilters2D(filtersize, maxscale=8, frameFunction=held, type=Float64)
    # Here, just normal dyadic dilations are considered and no higher forms. Also, the Filters in the final vector are stored from lowest to highest scale
    # The lowpass and the highpass filter are continious to the edge
    # good values for filters of real images which have been found by experiments. These K are scaling exponents which are used to scale the radius component. 

    # in case of non-dyadic filtering (not yet considered in julia)
    # possible values for a detailed decomposition (less and coarser is possible)
    # resolution, r, k
    # 32x32, 1, 5:-1:-1 
    # 64x64, 2, 11:-1:0 
    # 128x128, 3, 17:-1:1
    # 256x256, 4, 23:-1:2
    # 512x512, 5, 29:-1:3

    K = type.(maxscale-2:-1:-1)
    # only the pure radius is needed
    radius, _ = generateRadius2D(filtersize, type)
    filterbank = Array{Array{type}}(undef, length(K))
    #TODO benchmark this
    # the next loop is optimized to run also in parallel, but in case of small filterbanks, this is mostly not necessary and causes higher computation 
    # times through communications
    for scale = 1:length(K)
        #Threads.@threads for scale = 1:length(K)
        filterbank[scale] = frameFunction.(2^K[scale] .* radius)
    end
    filterbank[1] = (radius .* 2^K[1] .<= 0.25) + frameFunction.(2^K[1] .* radius) .* (radius .* 2^K[1] .> 0.25) # lowpass
    filterbank[length(K)] = 1 * (radius .* 2^K[end] .>= 0.25) + frameFunction.(2^K[end] .* radius) .* (radius .* 2^K[end] .< 0.25) # highpass
    return filterbank
end

"""
GenerateRieszFilters

Generate the Riesz filterbank for the according filtersize and depth of Riesz
"""
function generateRieszFilters2D(filtersize, maxrieszorder, type=Float64)
    # fast track for the simple case. This happens for ordinary Riesz Applications
    if maxrieszorder == 1
        filterbank = Array{Array{Complex{type}}}(undef, 2)
        radius, mesh_x, mesh_y, center_x, center_y = generateRadius2D(filtersize, type)
        filterbank[1] = mesh_x .* 1im ./ radius
        filterbank[1] = mesh_y .* 1im ./ radius
    end
    # calculate the size of the filter pyramid
    # it is guaranteed that this will be an Int
    rieszfilteramount = Int(0.5 * (maxrieszorder + 1) * (maxrieszorder + 2) - 1)
    # malloc arrays and create helper vars
    filterbank = Array{Array{Complex{type}}}(undef, rieszfilteramount)
    radius, mesh_x, mesh_y, center_x, center_y = generateRadius2D(filtersize, type)
    # create filterbank filter by filter and save it into the filterbank
    # the indexing method seems not consistend, but it fits together
    # Another approach would be to use  complex indexing inside a Pascal triangle, which is slower (benchmarked already)
    rieszfilteridx = 0
    for rieszorder = 1:maxrieszorder
        #complexify radius and prepare for division to have just one division per loop. Maybe more optimization is possible here
        radius_factor = ((1 ./ radius) .* (-1im)) .^ rieszorder
        # Riesz indices on this pyramid level
        riesz_idx_x = rieszorder:-1:0
        riesz_idx_y = rieszorder .- riesz_idx_x
        # level starts
        for generation_idx = 1:rieszorder+1
            rieszfilteridx = rieszfilteridx + 1
            filterbank[rieszfilteridx] = ((mesh_x .^ riesz_idx_x[generation_idx]) .* (mesh_y .^ riesz_idx_y[generation_idx])) .* radius_factor
            # deal with the center pixel which must be 0 by limit process.
            filterbank[rieszfilteridx][center_x, center_y] = zero(type)
        end
    end
    return filterbank
end

"""
GenerateRadius2D(filtersize)

Generates the Radius for a given 2D filtersize
"""

function generateRadius2D(filtersize, type=Float64)
    # malloc
    radius_x = generateRadius1D(filtersize[1], type)
    radius_y = generateRadius1D(filtersize[2], type)
    mesh_x = radius_x' .* ones(filtersize)
    mesh_y = radius_y .* ones(filtersize)
    radius = abs.(mesh_x + im * mesh_y)
    center_x = Int(floor(filtersize[1] / 2 + 1))
    center_y = Int(floor(filtersize[2] / 2 + 1))
    return (radius, mesh_x, mesh_y, center_x, center_y)
end #GenerateRadius2D

"""
GenerateRadius1D(filtersize)

Generates the Radius for a given 1D filtersize
"""
function generateRadius1D(filtersize, type=Float64)
    radius = range(-0.5, 0.5, filtersize + 1)
    radius = radius[1:end-1]
    return type.(radius)
end #GenerateRadius1D
