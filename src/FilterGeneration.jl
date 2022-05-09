"""
  GenerateRadius(filtersize)

Generates the Radius for a given filtersize
"""
function GenerateRadius2D(filtersize, type=Float64)
  # malloc
  finalRadius = Array{type}(undef, filtersize)
  radius_x = range(-0.5, 0.5, filtersize[1] + 1)
  radius_x = radius_x[1:end-1]
  radius_y = range(-0.5, 0.5, filtersize[2] + 1)
  radius_y = radius_y[1:end-1]

  mesh_x = radius_x' .* ones(filtersize)
  mesh_y = ones(filtersize) .* radius_y
  radius = abs(mesh_x + im * mesh_y)
  center_x = floor(filtersize(1) / 2 + 1)
  center_y = floor(filtersize(2) / 2 + 1)

  return (radius, mesh_x, mesh_y, center_x, center_y)
end #GenerateRadius2D

function GenerateRadius1D(filtersize, type=Float64)

end #GenerateRadius1D