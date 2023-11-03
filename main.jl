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
Packages = ["Profile", "CSV", "DataFrames", "Plots", "DataFrames", "Statistics", "Test", "LaTeXStrings"]
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

include("src\\integrators.jl")
include("src\\tools\\data_parser.jl")
include("src\\tools\\other_tools.jl")
include("test\\tests.jl")

run_tests()

function simulate_system(bodys::Vector{Vector{Any}}, force::Function, integrand=integ_euler, tot_time=50, η=0.01; dyn_dt::Bool=true, frames=100)
  E_tot_0, Ep_tot_0, Ev_tot_0 = tot_energy(bodys)
  P_tot_0 = tot_momentum(bodys)
  L_tot_0 = tot_ang_mom(bodys)
  t = 0
  
  # It is not yet clear, how many iterations we will have, so the array size is not predetermined
  E_n  = [E_tot_0]
  Ep_n = [Ep_tot_0]
  Ev_n = [Ev_tot_0]
  P_n  = [P_tot_0]
  L_n  = [L_tot_0]
  t_n  = [0.0]
  history = []
  history_t = []

  println("simulation started")
  progress = -1

  # Begin the simulation process
  @time "simulation time: " while t < tot_time
    # every 10 procent, display the progress:
    if t/tot_time > progress + 0.1
      progress = floor(t/tot_time, digits=1)
      println("Progress: ", progress * 100, "%")
    end

    # save the current state of the system, for future animation
    if size(history, 1) / frames <= t / tot_time
      append!(history, [deepcopy(bodys)])
      append!(history_t, [t])
    end
    
    # calulate time step for current iteration
    if !dyn_dt
      dt = η
    else
      dt = time_step(bodys, η)
    end

    bodys_copy = deepcopy(bodys)
    # update the bodys
    for i in axes(bodys, 1)
      f1(body) = force(body, bodys_copy[vcat(1:i-1, i+1:end)])
      # integrand(bodys[i], f1, dt)
      bodys[i] = integrand(bodys[i], f1, dt)
    end
    t += dt

    
    # Save current macroscopic characteristics
    En = tot_energy(bodys)
    append!(t_n, [t])
    append!(E_n, [En[1]])
    append!(Ep_n, [En[2]])
    append!(Ev_n, [En[3]])
    append!(P_n, [tot_momentum(bodys)])
    append!(L_n, [tot_ang_mom(bodys)])
  end

  println("simulation complete after $(size(t_n, 1)-1) iterations")

  metrik(vec) = collect(vec) .^ 2 |> sum |> sqrt
  P_abs = norm.(P_n)
  L_abs = norm.(L_n)

  function visualise_space(bodys, scope=-1; title="time = $(history_t[i])")
    # display(bodys)
    x = [b[1][1] for b in bodys]
    y = [b[1][2] for b in bodys]
    # z = [b[1][3] for b in bodys] |> ustrip
    plotsize = (750, 750)
    if scope==-1
      plot1 = scatter(x, y, title=title, size=plotsize, dpi=150)
    else
      plot1 = scatter(x, y, title=title, xlim=(scope[1][1], scope[1][2]), ylim=(scope[2][1], scope[2][2]), size=plotsize, dpi=150)
    end
    return plot1
  end

  
  
  plot1 = plot(t_n, [E_n, Ep_n, Ev_n], label=["E" "E_pot" "E_vel"], yscale=:identity, title="Energy", xlabel="Time", ylabel="Energy", dpi=300)
  # plot1 = plot(t_n, E_n, label="E", yscale=:identity, title="Energy", xlabel="Time", ylabel="Energy", dpi=300)
  plot2 = plot(t_n, P_abs, yscale=:identity, title="Momentum", xlabel="Time", ylabel="Momentum", dpi=300)
  plot3 = plot(t_n, L_abs, yscale=:identity, title="Anglular Momentum", xlabel="Time", ylabel="Anglular Momentum", dpi=300)
  plt   = plot(plot1, plot2, plot3, layout=(2,2), dpi=300)
  savefig(plt, "media\\plot")
  
  anim2 = @animate for i in eachindex(history)
    gridsize = 5
    time = history_t[i] |> string
    time = time * "0000000"
    time = time[1:(div(length(time), 10) + 3)]
    visualise_space(history[i], ((-gridsize, gridsize), (-gridsize, gridsize)), title="time = $(time)")
  end
  
  path = string(@__DIR__) * "\\media\\GIF.gif"
  println("Saving gif")
  gif(anim2, path)
end

function evaluate_2body_problem()
  # define force (gravity)
  function gravity(body, bodys)
    G = 1
    F = (0, 0, 0)
    r_min = 0.005

    for b in bodys
      r = norm(b[1] .- body[1])
      if r > r_min
        F = F .- G * b[3] * body[3] / r^3 .* (body[1] .- b[1])
      else
        F = F .- G * b[3] * body[3] / (r_min^2 * r) .* (body[1] .- b[1])
      end
    end

    return F
  end

  path_2body = "data\\2_body.csv"
  bodys_original = parse_data(path_2body)
  center_system!(bodys_original)
  norm_mass!(bodys_original)
  
  # Parameters for each integration method:
  integrators   = [integ_euler, integ_euler_cromer, integ_rk2, integ_rk4]
  time_scale    = 2pi .* [200, 200, 200, 200] .* 0.25
  stepsize_mult = [1, 1, 1, 1]
  dynamic_step  = [true, true, true, true]
  plotname      = ["2-Body_Euler", "2-Body_Euler-Cromer", "2-Body_RK2", "2-Body_RK4"]
  integ_name    = ["Euler", "Euler-Cromer", "RK2", "RK4"]
  # containers for Plots
  plot_E = plot(xlabel="Orbits", title=L"$\log_{10}|\frac{E-E_0}{E_0}|$", dpi=300)
  plot_j = plot(xlabel="Orbits", title=L"|\vec{j}|", dpi=300)
  plot_e = plot(xlabel="Orbits", title=L"\log_{10}|\frac{e-e_0}{e_0}|", dpi=300)
  plot_a = plot(xlabel="Orbits", title=L"\log_{10}|\frac{a_e-a_{e,0}}{a_{e,0}}|", dpi=300)

  for i in integrators |> eachindex
    bodys = deepcopy(bodys_original) # deepcopy bodys each time, to avoid aliasing

    println("\nEvaluation for integrator ", plotname[i], " started")
    t = 0.0
    t_end = time_scale[i]
    η = stepsize_mult[i] * 10e-3
    integrator = integrators[i]

    # containers for the macroscopic data we want to track
    E_n   = [tot_energy(bodys)[1]]
    j_n   = [spec_ang_mom(bodys[1])]
    e_n   = [runge_lenz_vec(bodys[1])]
    a_e_n = [large_semi_axis(bodys[1])]
    t_n   = [t]


    progress = -1

    while t <= t_end
      # output progress every 10%
      if t/t_end > progress + 0.1
        progress = floor(t/t_end, digits=1)
        println("Progress: ", round(Int, progress * 100), "%")
      end

      dt = dynamic_step[i] ? time_step(bodys, η) : η # calculate current time step, dependant on curvature of particles
      t += dt
      # update each body
      start_bodys = deepcopy(bodys)
      for i in eachindex(bodys)
        f1(body) = gravity(body, start_bodys[vcat(1:i-1, i+1:end)])
        bodys[i] = integrator(bodys[i], f1, dt)
      end
      append!(E_n, [tot_energy(bodys)[1]])
      append!(j_n, [spec_ang_mom(bodys[1])])
      append!(e_n, [runge_lenz_vec(bodys[1])])
      append!(a_e_n, [large_semi_axis(bodys[1])])
      append!(t_n, [t])
    end
    
    # Plot the makroscopic data:
    # y1 = E_n
    # y2 = norm.(j_n)
    # y3 = norm.(e_n)
    # y4 = a_e_n

    y1 = log10.(abs.((E_n .- E_n[1]) ./ E_n[1]))
    y2 = norm.(j_n)
    y3 = log10.(abs.((norm.(e_n) .- norm.(e_n)[1])/norm.(e_n)[1]))
    y4 = log10.(abs.((a_e_n .- a_e_n[1])/a_e_n[1]))

    linestyle = :solid
    if integ_name[i] == "RK4"
      linestyle = :dashdotdot
    end

    plot_E = plot!(plot_E, t_n ./ 2pi, y1, label=integ_name[i], linestyle=linestyle)
    plot_j = plot!(plot_j, t_n ./ 2pi, y2, label="", linestyle=linestyle)
    plot_e = plot!(plot_e, t_n ./ 2pi, y3, label="", linestyle=linestyle)
    plot_a = plot!(plot_a, t_n ./ 2pi, y4, label="", linestyle=linestyle)
  end

  println("generating Plots...")
  tot_plot = plot(plot_E, plot_j, plot_e, plot_a, layout=(2,2), dpi=300, size=(600, 400) .* 1.3)
  savefig(tot_plot, "media\\total_plot")
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
    r_min = 0.005
    if size(bodys, 1) > 20
      Threads.@threads for b in bodys
        @fastmath r = norm(b[1] .- body[1])
        if r > r_min
          @fastmath F = F .- G * b[3] * body[3] / r^3 .* (body[1] .- b[1])
        else
          @fastmath F = F .- G * b[3] * body[3] / (r_min^2 * r) .* (body[1] .- b[1])
        end
      end
    else
      for b in bodys
        @fastmath r = norm(b[1] .- body[1])
        if r > r_min
          @fastmath F = F .- G * b[3] * body[3] / r^3 .* (body[1] .- b[1])
        else
          @fastmath F = F .- G * b[3] * body[3] / (r_min^2 * r) .* (body[1] .- b[1])
        end
      end
    end
    return F
  end

  # simulate_system(data1, gravity, integ_euler, 25, 0.001, dyn_dt=false)
  simulate_system(data1, gravity, integ_rk4, 2pi*200, 0.001, dyn_dt=false, frames=250)

  return nothing
end

if string(@__MODULE__) == "Main"
  # start()
  evaluate_2body_problem()
end