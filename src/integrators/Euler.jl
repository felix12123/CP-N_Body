# Es folgt nun die Implementation der Integratoren





# ****************************
# ***** Euler-Integrator *****
# ****************************

function integ_euler(body, force, dt)

	# Berechne die neuen Orte und Geschwindigkeiten:
	x1 = body[1] .+ (body[2] .* dt)							# x1 = x + v*dt
	v1 = body[1] .+ ((force(body) ./ body[3]) .* dt)		# v1 = v + a*dt = v + (f/m)*dt

	# Schreibe neue Werte in body:
	body[1] = x1
	bdoy[2] = v1

	# Gebe body zurück:
	return body
end





# *************************
# ***** Runge-Kutta 2 *****
# *************************

function integ_rk2(body, force, dt)
	# x1 = x0 + dt/2 * (v(t0) + v(t0 + dt))
	# v1 = v0 + dt/2 * (a(x) + a(x0 + dx))

	# Berechne x1, x2, v1, v2 für die spätere Mittelung.
	# Hier nun x1 und v1:
	x1 = body[1] .+ (body[2] .* dt)
	v1 = body[1] .+ ((force(body) ./ body[3]) .* dt)

	# Schreibe neue Werte in body, damit neue 
	# Beschleunigung berechnet werden kann:
	body[1] = x1
	bdoy[2] = v1

	# Hier nun x2 und v2:
	x2 = body[1] .+ (body[2] .* dt)
	v2 = body[1] .+ ((force(body) ./ body[3]) .* dt)

	# Führe nun Mittelung aus und schreibe diese
	# Werte in den finalen body:
	xnew = mean.(x1, x2)
	vnew = mean.(v1, v2)

	# Schreibe gemittelte Werte in body rein
	# für die finale Ausgabe:
	body[1] = xnew
	bdoy[2] = vnew

	return body
end





# *************************
# ***** Runge-Kutta 2 *****
# *************************