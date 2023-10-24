
# Normate Mass of system
function norm_mass!(bodys::Vector{Vector{Any}})
  if !isa(bodys[1][3], Float64) @warn "norm_mass: Type of Data not valid for normation of mass" end

  M = 0
  for i in axes(bodys, 1)
    M += bodys[i][3]
  end
  
  for i in axes(bodys, 1)
    bodys[i][3] /= M
  end
end