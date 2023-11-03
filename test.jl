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


function fast_sinus()
  for i in 1:10000000
    @fastmath(sin(log(i)))
  end
end
function sinus()
  for i in 1:10000000
    sin(log(i))
  end
end

function start()
  
  @time "fast" fast_sinus()
  @time "fast" fast_sinus()
  @time "not fast" sinus()
  @time "not fast" sinus()

end
start()