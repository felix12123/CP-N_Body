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

include("src\\integrators\\Euler.jl")
include("src\\tools\\data_parser.jl")
include("src\\tools\\other_tools.jl")
include("test\\integrator_tests.jl")

run_tests()



function simulate_system(bodys::Vector{Vector{Any}}, force::Function, integrand=integ_euler, tot_time=100, η=0.05, )
  E_tot = tot_energy(bodys)
  P_tot = tot_momentum(bodys)
  t = 0
  E_t = []
  while t < tot_time
    dt = time_step(bodys, η)
    t += dt
    for i in axis(bodys, 1)
      acc(pos) = force(pos, bodys[vcat(1:i-1, i+1:end)]) / bodys[i][3]
      bodys[i][1] = bodys[i][2] .* dt
      bodys[i][2] = integrand(acc, bodys[i][2], dt)
    end
  end
end


function start()
  # Read starting parameters from data files
  data1 = parse_data("data\\2_body.csv")
  data2 = parse_data("data\\3_body.csv")
  data3 = parse_data("data\\100_body.csv")
  data4 = parse_data("data\\1000_body.csv")

  # simulate_system(data1, 10000)

  return nothing
end


if string(@__MODULE__) == "Main"
  @time "whole progarm duration" start()
end

