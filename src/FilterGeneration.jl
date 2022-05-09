"""
  GenerateRadius(filtersize)

Generates the Radius for a given filtersize of any size
"""
function GenerateRadius
  
end #GenerateRadius

"""
  GenerateRadius2D(filtersize)

Generates the Radius for a given 2D filtersize
"""
function GenerateRadius2D(filtersize, type=Float64)
  # malloc
  radius_x = GenerateRadius1D(filtersize[1], type)
  radius_y = GenerateRadius1D(filtersize[2], type)
  mesh_x = radius_x' .* ones(filtersize)
  mesh_y = radius_y .* ones(filtersize)
  radius = abs.(mesh_x + im * mesh_y)
  center_x = floor(filtersize[1] / 2 + 1)
  center_y = floor(filtersize[2] / 2 + 1)
  return (radius, mesh_x, mesh_y, center_x, center_y)
end #GenerateRadius2D

"""
  GenerateRadius1D(filtersize)

Generates the Radius for a given 1D filtersize
"""
function GenerateRadius1D(filtersize, type=Float64)
  radius = range(-0.5, 0.5, filtersize + 1)
  radius = radius[1:end-1]
  return type.(radius)
end #GenerateRadius1D