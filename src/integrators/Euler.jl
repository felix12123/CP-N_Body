# Es folgt nun die Implementation der Integratoren





# ****************************
# ***** Euler-Integrator *****
# ****************************

function integ_euler(body0, force, dt)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Berechne die neuen Orte und Geschwindigkeiten:
	x1 = body[1] .+ (body[2] .* dt)								# x1 = x + v*dt
	v1 = body[1] .+ ((force(body) ./ body[3]) .* dt)			# v1 = v + a*dt = v + (f/m)*dt

	# Schreibe neue Werte in body:
	body[1] = x1
	body[2] = v1

	# Gebe body zurück:
	return body
end

# Euler-Integrator, bei dem die Eingabe überschrieben wird (benutzen wir wahrschienlich nicht):
function integ_euler!(body, force, dt)

	# Berechne die neuen Orte und Geschwindigkeiten:
	x1 = body[1] .+ (body[2] .* dt)							# x1 = x + v*dt
	v1 = body[1] .+ ((force(body) ./ body[3]) .* dt)		# v1 = v + a*dt = v + (f/m)*dt

	# Schreibe neue Werte in body:
	body[1] = x1
	body[2] = v1
end





# *************************
# ***** Runge-Kutta 2 *****
# *************************

function integ_rk2(body0, force, dt)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Speichere x0, v0 und a0:
	x0 = body[1]
	v0 = body[2]
	a0 = (force(body) ./ body[3])

	# -Berechne die benötigten Werte:
	x1 = body[1] .+ (body[2] .* dt)
	
	# Schreibe neue Werte in body, damit neue Beschleunigung berechnet werden kann:
	body[1] = x1

	# Weiter in der Berechnung:
	a1 = (force(body) ./ body[3])
	v1 = body[2] .+ (a1 .* dt)

	# Berechne Vorhersage nach RK2:
	xnew = x0 + (dt/2) * mean.(v0, v1)
	vnew = v0 + (dt/2) * mean.(a0, a1)

	# Schreibe gemittelte Werte in body rein für das finale Ergebnis:
	body[1] = xnew
	body[2] = vnew
end





# *************************
# ***** Runge-Kutta 4 *****
# *************************

function integ_rk4(body0, force, dt)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Häufig verwendete Elemente sollten nicht durchweg auf's neue berechnet werden müssen:
	m = body[3]
	dt_2 = dt/2

	# Rechnungen zu K1:
	x0 = body[1]										# Ziehe Anfangswert von x
	v0 = body[2]										# Ziehe Anfangswert von v
	a0 = (force(body) ./ m)								# Beschleunigung am Startpunkt

	# Rechnungen zu K2
	x1 = x0 .+ (v0 .* dt_2)								# Berechne x1 via: x1 = x0 + v0*(dt/2)
	v1 = v0 .+ (a0 .* dt_2)								# Berechne v1 via: v1 = v0 + a0*(dt/2)
	body[1] = x1										# Schreibe x1 in body, damit a1 berechnet werden kann
	body[2] = v1										# Schreibe v1 in body für spätere Berechnungen von x2
	a1  = (force(body) ./ m)							# Berechne a1 über die Kraft: f/m = a

	# Ermittel nun K3:
	x2 = x0 .+ (v1 .* dt_2)								# x2 = x0 + v1*(dt/2)
	v2 = v0 .+ (a1 .* dt_2)								# Berechne v2 via: v2 = v0 + a1*(dt/2)
	body[1] = x2										# Schreibe x2 in body, damit a2 berechnet werden kann
	body[2] = v2										# Schreibe v2 in body für spätere Berechnungen von x3
	a2 = (force(body) ./ m)								# Berechne a2 über die Kraft: f/m = a

	# Ermittel nun K4:
	x3 = x0 .+ (v2 .* dt)								# Berechne x3 via: x3 = x0 + v2*dt
	v3 = v0 .+ (a2 .* dt)								# Berechne v3 via: v3 = v0 + a2*dt
	body[1] = x3										# Schreibe x3 in body, damit a3 berechnet werden kann
	body[2] = v3										# Schreibe v3 in body
	a3 = (force(body) ./ m)								# Berechne a3 über die Kraft: f/m = a

	# Berechne nun Vorhersage nach RK4:
	xnew = x0 + (dt/6)*mean.(v0, v1, v2, v3)
	vnew = v0 + (dt/6)*mean.(a0, a1, a2, a3)

	# Schreibe gemittelte Werte in body rein für das finale Ergebnis:
	body[1] = xnew
	body[2] = vnew

	# Gebe Ergebnis-body aus:
	return body
end