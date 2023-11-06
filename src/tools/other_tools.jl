using LinearAlgebra
import LinearAlgebra.cross
cross(a::NTuple{3, Float64}, b::NTuple{3, Float64}) = Tuple(cross(collect(a), collect(b)))


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
  abs_acc = Vector{Float64}(undef, size(bodys, 1)) # Vector to store the acceleration of each body
  abs_acc_change = Vector{Float64}(undef, size(bodys, 1)) # Vector to store the change in acceleration of each body
  
  for i in eachindex(abs_acc)
    a_i = (0, 0, 0)
    a_dot_i = (0, 0, 0)

    for j in eachindex(abs_acc)
      if i != j
        r_ij = bodys[i][1] .- bodys[j][1]
        abs_r_ij = norm(r_ij)
        v_ij = bodys[i][2] .- bodys[j][2]
        
        a_i = a_i .+ bodys[j][3] / abs_r_ij^3 .* r_ij

        a_dot_i = a_dot_i .+ bodys[j][3] .* (v_ij ./ abs_r_ij^3 .- 3/abs_r_ij^5 * dot(v_ij, r_ij) .* r_ij)
      end
    end
    abs_acc[i] = norm(a_i)
    abs_acc_change[i] = norm(a_dot_i)
  end

  step = dt0 * minimum(abs_acc ./ abs_acc_change)
  if step == 0 @error "Time step evaluated to 0" end
  return min(step, 15*dt0)
end

# calculates the total energy of the system
function tot_energy(bodys::Vector{Vector{Any}})
  G = 1
  E_pot_tot = 0
  E_vel_tot = 0
  for b in bodys
    E_vel = 1/2 * b[3] * sum(b[2] .^ 2)
    E_pot = 0
    for b1 in bodys
      if b1 != b
        E_pot += G * b1[3] * b[3] / sqrt(sum((b1[1] .- b[1]) .^ 2))
      end
    end
    E_vel_tot += E_vel
    E_pot_tot += E_pot
  end
  return E_vel_tot - E_pot_tot, E_pot_tot, E_vel_tot
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
    L = L .+ cross(b[1], b[2]) .* b[3]
  end

  return L
end

function spec_ang_mom(body)
  body[2] × body[1]
end

function runge_lenz_vec(body)
  body[2] × spec_ang_mom(body) .+ body[1] ./ norm(body[1])
end


function large_semi_axis(body)
  j = spec_ang_mom(body)
  e = runge_lenz_vec(body)
  return dot(j, j) * (1-dot(e, e))
end