"""
    GenerateRieszFilters(filtersize[, maxscale, framefunction, type])

Generate the Frame filterbank for the according filtersize and choosen Motherwavelet

# Arguments
- `maxscale`: blabla
- `framefunction`: bla
- `type`: some info about type
"""
function GenerateFrameFilters2D(filtersize, maxscale=8, framefunction=Held, type=Float64)
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

    decompositionscheme = type.(maxscale-2:-1:-1)
    # only the pure radius is needed
    radius, _ = GenerateRadius2D(filtersize, type)
    filterbank = Array{Array{type}}(undef, length(decompositionscheme))
    #TODO benchmark this
    # the next loop is optimized to run also in parallel, but in case of small filterbanks, this is mostly not necessary and causes higher computation 
    # times through communications
    # leave the high- and lowpass filter for later computation (enlargement)
    for scale = 2:length(decompositionscheme)-1
        #Threads.@threads for scale = 1:length(K)
        filterbank[scale] = framefunction.(2^decompositionscheme[scale] .* radius)
    end
    # enlarge the high- and lowpass filter at the max point
    filterbank[1] =
        (radius .* 2^decompositionscheme[1] .<= 0.25) +
        framefunction.(2^decompositionscheme[1] .* radius) .* (radius .* 2^decompositionscheme[1] .> 0.25) # lowpass
    filterbank[length(decompositionscheme)] =
        1 * (radius .* 2^decompositionscheme[end] .>= 0.25) +
        framefunction.(2^decompositionscheme[end] .* radius) .* (radius .* 2^decompositionscheme[end] .< 0.25) # highpass
    return filterbank
end

"""
GenerateRieszFilters2D(filtersize[, maxrieszorder, type])

Generate the Riesz filterbank for the according filtersize and depth of Riesz
"""
function GenerateRieszFilters2D(filtersize, maxrieszorder=1, type=Float64)
    # fast track for the simple case. This happens for ordinary Riesz Applications

    if maxrieszorder == 1
        filterbank = Array{Array{Complex{type}}}(undef, 2)
        radius, mesh_x, mesh_y, center_x, center_y = GenerateRadius2D(filtersize, type)
        filterbank[1] = mesh_x .* -1im ./ radius
        filterbank[1][center_x, center_y] = zero(type)
        filterbank[2] = mesh_y .* -1im ./ radius
        filterbank[2][center_x, center_y] = zero(type)
        return filterbank
    end

    # This part is for higher order Riesz transforms only. Not useful in case of ordinary image processing

    # calculate the size of the filter pyramid
    # it is guaranteed that this will be an Int
    rieszfilteramount = Int(0.5 * (maxrieszorder + 1) * (maxrieszorder + 2) - 1)
    # malloc arrays and create helper vars
    filterbank = Array{Array{Complex{type}}}(undef, rieszfilteramount)
    radius, mesh_x, mesh_y, center_x, center_y = GenerateRadius2D(filtersize, type)
    # create filterbank filter by filter and save it into the filterbank
    # the indexing method seems not consistend, but it fits together
    # Another approach would be to use  complex indexing inside a Pascal triangle, which is slower (benchmarked already)
    rieszfilter_idx = 0
    for rieszorder = 1:maxrieszorder
        #complexify radius and prepare for division to have just one division per loop. Maybe more optimization is possible here
        radiusfactor = ((1 ./ radius) .* (-1im)) .^ rieszorder
        # Riesz indices on this pyramid level
        riesz_idx_x = rieszorder:-1:0
        riesz_idx_y = rieszorder .- riesz_idx_x
        # level starts
        for generation_idx = 1:rieszorder+1
            rieszfilter_idx = rieszfilter_idx + 1
            filterbank[rieszfilter_idx] = ((mesh_x .^ riesz_idx_x[generation_idx]) .* (mesh_y .^ riesz_idx_y[generation_idx])) .* radiusfactor
            # deal with the center pixel which must be 0 by limit process.
            filterbank[rieszfilter_idx][center_x, center_y] = zero(type)
        end
    end
    return filterbank
end

"""
GenerateRadius2D(filtersize[, type])

Generates the Radius for a given 2D filtersize
"""
function GenerateRadius2D(filtersize, type=Float64)
    # malloc
    radius_x = GenerateRadius1D(filtersize[1], type)
    radius_y = -GenerateRadius1D(filtersize[2], type) # because of orientation on the axis and on images
    mesh_x = radius_x' .* ones(filtersize)
    mesh_y = radius_y .* ones(filtersize)
    radius = abs.(mesh_x + im * mesh_y)
    center_x = Int(floor(filtersize[1] / 2 + 1))
    center_y = Int(floor(filtersize[2] / 2 + 1))
    return (radius, mesh_x, mesh_y, center_x, center_y)
end #GenerateRadius2D

"""
GenerateRadius1D(filtersize[, type])

Generates the Radius for a given 1D filtersize
"""
function GenerateRadius1D(filtersize, type=Float64)
    radius = range(-0.5, 0.5, filtersize + 1)
    radius = radius[1:end-1]
    return type.(radius)
end #GenerateRadius1D
