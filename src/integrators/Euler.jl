function integ_euler(v, x0, t0, dt)
	x1 = x0 + v(t0)*dt
end


function integ_rk2(v, x0, t0, dt)
	# x1 = x0 + dt/2 * (v(t0) + v(t0 + dt))
	
	# v1 = v0 + dt/2 * (a(x) + a(x0 + dx))
end

