using CSV, DataFrames


# Parses csv files of the format x;y;z;v_x;v_y;v_z;m into a Vector of Bodys
function parse_data(datafile::String; delim=";")
  data = CSV.File(datafile, delim=delim, header=false) |> DataFrame
  rename!(data, [:x, :y, :z, :vx, :vy, :vz, :m])
  
  len = size(data, 1)
  bodys = Vector{typeof([(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0, 0.0])}(undef, len)
  for i in 1:len
    bodys[i] = [(data.x[i], data.y[i], data.z[i]), (data.vx[i], data.vy[i], data.vz[i]), data.m[i], 0]
  end

  return bodys
end

