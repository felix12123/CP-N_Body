using LinearAlgebra


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

# transforms the system to the center of gravity system, so it doesnt drift
function center_system!(bodys::Vector{Vector{Any}})
  pos_correction = (0.0, 0.0, 0.0)
  vel_correction = (0.0, 0.0, 0.0)
  M = 0

  for i in axes(bodys, 1)
    pos_correction = pos_correction .- (bodys[i][1] .* bodys[i][3])
    vel_correction = vel_correction .- (bodys[i][2] .* bodys[i][3])
    M += bodys[i][3]
  end

  pos_correction = pos_correction ./ M
  vel_correction = vel_correction ./ M

  for i in axes(bodys, 1)
    bodys[i][1] = bodys[i][1] .+ pos_correction
    bodys[i][2] = bodys[i][2] .+ vel_correction
  end
  
  return nothing
end


# calculates the best time step size
# FIXME, Definition falsch verstanden, to be implemented
function time_step(acc1, acc2, dt0)
  return dt0
  abs_vel1 = [sqrt(sum(body[i][2] .^2)) for body in bodys1]
  abs_vel2 = [sqrt(sum(body[i][2] .^2)) for body in bodys2]
  # abs_pos = [sqrt(sum(body[i][1] .^2)) for body in bodys]
  # step = dt0 * minimum(abs_pos ./ abs_vel)
  # if step == 0 @error "Time step evaluated to 0" end
  # return step
end


# calculates the total energy of the system
function tot_energy(bodys::Vector{Vector{Any}})::Float64
  G = 1
  E = 0
  for b in bodys
    E_vel = 1/2 * b[3] * sum(b[2] .^ 2)
    E_pot = 0
    for b1 in bodys
      if b1 != b
        E_pot += G * b1[3] * b[3] / sqrt(sum((b1[1] .- b[1]) .^ 2))
      end
    end
    E += E_vel - E_pot
  end
  return E
end


# calculates the total impulse of the system
function tot_momentum(bodys::Vector{Vector{Any}})::Tuple{Float64, Float64, Float64}
  P = (0.0, 0.0, 0.0)
  
  for b in bodys
    P = P .+ b[3] .* b[2]
  end

  return P
end


# calculates the total angular momentum of the system
function tot_ang_mom(bodys::Vector{Vector{Any}})::Tuple{Float64, Float64, Float64}
  L = (0.0, 0.0, 0.0)
  
  for b in bodys
    L = L .+ Tuple(cross(collect(b[1]), collect(b[2])) .* b[3])
  end

  return L
end