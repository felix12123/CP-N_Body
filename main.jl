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
Packages = ["CSV", "DataFrames"]
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

include(string(@__DIR__) * "\\src\\integrators\\Euler.jl")
include(string(@__DIR__) * "\\src\\tools\\data_parser.jl")
include(string(@__DIR__) * "\\test\\integrator_tests.jl")

run_tests()



