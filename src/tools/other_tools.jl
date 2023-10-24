
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
function time_step(bodys, dt0)
  abs_pos = [sqrt(sum(body[i][1] .^2)) for body in bodys]
  abs_vel = [sqrt(sum(body[i][2] .^2)) for body in bodys]
  step = dt0 * minimum(abs_pos ./ abs_vel)
  if step == 0 @error "Time step evaluated to 0" end
  return step
end


# calculates the total energy of the system
function tot_energy(bodys::Vector{Vector{Any}})::Float64
  E = 0
  for b in bodys
    E += 1/2 * b[3] * sum(b[2] .^ 2)
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