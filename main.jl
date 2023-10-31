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
Packages = ["CSV", "DataFrames", "Plots", "DataFrames", "Statistics", "Test"]
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


function simulate_system(bodys::Vector{Vector{Any}}, force::Function, integrand=integ_euler, tot_time=50, η=0.01; dyn_dt::Bool=true)
  E_tot_0 = tot_energy(bodys)
  P_tot_0 = tot_momentum(bodys)
  L_tot_0 = tot_ang_mom(bodys)
  t = 0
  
  # It is not yet clear, how many iterations we will have, so the array size is not predetermined
  E_n = [E_tot_0]
  P_n = [P_tot_0]
  L_n = [L_tot_0]
  t_n = [0.0]
  history = []
  # We need to save the previous accelerations of our system, to calculate the dynamic time step
  # prev_accs = undef
  # prev_vels = undef
  # accs = undef

  println("simulation started")
  # Begin the simulation process
  @time "simulation time: " while t < tot_time
    if size(t, 1) % 10 == 1
      append!(history, [deepcopy(bodys)])
    end
    # calulate time step for current iteration
    if t == 0 | !dyn_dt
      dt = η
    else
      # accs = 
      dt = η #time_step(prev_accs, accs, η)
    end

    # prev_vels = copy(bodys)
    
    # update the bodys
    for i in axes(bodys, 1)
      f1(body) = force(body, bodys[vcat(1:i-1, i+1:end)])
      bodys[i] = integrand(bodys[i], f1, dt)
    end
    t += dt

    # vels = [b[2] for b in bodys]
    # prev_accs =  (vels .- prev_vels) ./ dt
    
    # Save current macroscopic characteristics
    append!(t_n, [t])
    append!(E_n, [tot_energy(bodys)])
    append!(P_n, [tot_momentum(bodys)])
    append!(L_n, [tot_ang_mom(bodys)])
  end

  println("simulation complete after $(size(t_n, 1)) iterations")

  metrik(vec) = collect(vec) .^ 2 |> sum |> sqrt
  P_abs = metrik.(P_n)
  L_abs = metrik.(L_n)

  function visualise_space(bodys, scope=-1)
    # display(bodys)
    x = [b[1][1] for b in bodys]
    y = [b[1][2] for b in bodys]
    # z = [b[1][3] for b in bodys] |> ustrip
    plotsize = (750, 750)
    if scope==-1
      plot1 = scatter(x, y, size=plotsize, dpi=150)
    else
      plot1 = scatter(x, y, xlim=(scope[1][1], scope[1][2]), ylim=(scope[2][1], scope[2][2]), size=plotsize, dpi=150)
    end
    return plot1
  end

  
  
  plot1 = plot(t_n, E_n, yscale=:identity, title="Energy", xlabel="Time", ylabel="Energy", dpi=300)
  plot2 = plot(t_n, P_abs, yscale=:identity, title="Momentum", xlabel="Time", ylabel="Momentum", dpi=300)
  plot3 = plot(t_n, L_abs, yscale=:identity, title="Anglular Momentum", xlabel="Time", ylabel="Anglular Momentum", dpi=300)
  plt   = plot(plot1, plot2, plot3, layout=(2,2), dpi=300)
  display(plt)
  
  anim2 = @animate for i in eachindex(history)[1:div(end, 50):end]
    gridsize = 5
    visualise_space(history[i], ((-gridsize, gridsize), (-gridsize, gridsize)))
  end

  path = string(@__DIR__) * "\\GIF.gif"
  println("Saving gif in ", "TestGif")
  gif1 = gif(anim2, path)
end


function start()
  # Read starting parameters from data files
  data1 = parse_data("data\\2_body.csv")
  data2 = parse_data("data\\3_body.csv")
  data3 = parse_data("data\\100_body.csv")
  data4 = parse_data("data\\1000_body.csv")

  center_system!(data1)
  center_system!(data2)
  center_system!(data3)
  center_system!(data4)

  norm_mass!(data1)
  norm_mass!(data2)
  norm_mass!(data3)
  norm_mass!(data4)

  function gravity(body, bodys)
    G = 1
    F = (0, 0, 0)
    for b in bodys
      r = norm(b[1] .- body[1])
      if r > 0.001
        F = F .- G * b[3] * body[3] / r^3 .* (body[1] .- b[1])
      else

        F = F .- G * b[3] * body[3] / 0.001^2 .* (body[1] .- b[1]) ./ r
      end
        
    end
    return F
  end

  # simulate_system(data1, gravity, integ_euler, 25, 0.001, dyn_dt=false)
  simulate_system(data4, gravity, integ_euler, 5, 0.1, dyn_dt=false)

  return nothing
end


if string(@__MODULE__) == "Main"
  @time "whole program duration" start()
end

