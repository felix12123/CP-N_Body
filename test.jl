using Pkg
function installed()
  deps = Pkg.dependencies()
  installs = Dict{String, VersionNumber}()
  for (uuid, dep) in deps
    dep.is_direct_dep || continue
    dep.version === nothing && continue
    installs[dep.name] = dep.version
  end
  return installs
end
# Check if packages are installed, else install them
Packages = []
installed_Packages = keys(installed())
for Package in Packages
  if !(Package in installed_Packages)
    try
      eval(Meta.parse("using $Package"))
    catch
      println("Package $Package was not found. Installation started")
      Pkg.add(Package)
      eval(Meta.parse("using $Package"))
    end
  else
    eval(Meta.parse("using $Package"))
  end
end


function count_julia_lines(path::String)::Int
  if isfile(path)
    if split(path, ".")[end] == "au3"
      return open(path) do file
        lines = size(eachline(file) |> collect, 1)
        return lines
      end
      return 0
    end
  elseif isdir(path)
    return sum(count_julia_lines.(readdir(path, join=true)))
  end
  return 0
end

