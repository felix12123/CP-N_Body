# Es folgt nun die Implementation der Integratoren





# *****************
# ***** Euler *****
# *****************

function integ_euler(body0, force, dt)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Damit body[2] nicht 2x ausgelesen werden muss:
	v0 = body[2]

	# Berechne die neuen Orte und Geschwindigkeiten:
	v1 = v0 .+ ((force(body) ./ body[3]) .* dt)				# v1 = v + a*dt = v + (f/m)*dt
	x1 = body[1] .+ (v0 .* dt)								# x1 = x + v*dt

	# Schreibe neue Werte in body:
	body[1] = x1
	body[2] = v1

	# Gebe body zurück:
	return body
end

# Schnellere in situ Variante:
function fast_integ_euler!(body, force, dt)
	body[2] = body[2] .+ ((force(body) ./ body[3]) .* dt)
	body[1] = body[1] .+ (body[2] .* dt)
end





# ************************
# ***** Euler-Cromer *****
# ************************

function integ_euler_cromer(body0, force, dt)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Damit body[2] nicht 2x ausgelesen werden muss:
	v0 = body[2]

	# Berechne die neuen Orte und Geschwindigkeiten:
	v1 = v0 .+ ((force(body) ./ body[3]) .* dt)				# v1 = v0 + a0*dt = v + (f/m)*dt
	x1 = body[1] .+ ((1/2) .* (v0 .+ v1) .* dt)				# x1 = x0 + (1/2)*(v0 + v1)dt

	# Schreibe neue Werte in body:
	body[1] = x1
	body[2] = v1

	# Gebe body zurück:
	return body
end

# Schnellere in situ Variante:
function fast_integ_euler_cromer!(body, force, dt)
	v0 = body[2]
	body[2] = body[2] .+ ((force(body) ./ body[3]) .* dt)
	body[1] = body[1] .+ ((1/2) .* (v0 .+ body[2]) .* dt)
end





# ***************************
# ***** Velocity-Verlet *****
# ***************************

function integ_velocity_verlet(body0, force, dt)
	
	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Erforderliche Größen:
	# Speichere x0, v0 und a0:
	x0 = body[1]
	v0 = body[2]
	m =  body[3]
	a0 = (force(body) ./ m)

	# **********************************************************************
	# Velocity-Verlet für konstante Zeitschritte:
	#	# Berechnungsvorschrift für den neuen Ort:
	#	x1 = x0 .+ (v0 .* dt) + ((1/2) * a0 * dt^2)

	#	# Schreibe neues x (= x1) in body, um neues a (= a1) zu berechnen.
	#	# Dieses benötigt man um neues v (= v1) zu berechnen:
	#	body[1] = x1
	#	a1  = (force(body) ./ m)
	#	v1 = v0 .+ ((1/2) * mean.(a0, a1) * dt)

	#	# Schreibe neues v (= v1) in body:
	#	body[2] = v1

	#	# Gebe Ergebnis zurück:
	#	return body
	# **********************************************************************

	# Velocity-Verlet für variable Zeitschritte (kick-drift-kick-Methode):
	v1_05 = v0 .+ ((1/2) .* a0 .* dt)					# (1.38)
	x1 	  = x0 .+ (v1_05 .* dt)							# (1.39)

	#	# Schreibe neues x (= x1) in body, um neues a (= a1) zu berechnen.
	#	# Dieses benötigt man um neues v (= v1) zu berechnen:
	body[1] = x1
	a1  	= (force(body) ./ m)
	v1 		= v1_05 .+ (1/2 .* a1 .* dt)				# (1.40)

	# Schreibe neues v (= v1) in body:
	body[2] = v1

	# Gebe Ergebnis zurück:
	return body
end

# Schnellere (Kick-Drift-Kick-Variante) in situ Variante:
function fast_integ_velocity_verlet!(body, force, dt)
	v1_05 = body[2] .+ ((1/2) .* (force(body) ./ body[3]) .* dt)
	body[1] = body[1] .+ (v1_05 .* dt)	
	body[2]	= v1_05 .+ (1/2 .* (force(body) ./ body[3]) .* dt)
end






# *************************
# ***** Runge-Kutta 2 *****
# *************************

# Hier der RK2-Integrator bzw. Heun-Integrator

function integ_rk2(body0, force, dt)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Speichere x0, v0 und a0:
	x0 = body[1]
	v0 = body[2]
	m =  body[3]
	a0 = (force(body) ./ m)

	# -Berechne die benötigten Werte:
	x1 = x0 .+ (v0 .* dt)
	
	# Schreibe neue Werte in body, damit neue Beschleunigung berechnet werden kann:
	body[1] = x1

	# Weiter in der Berechnung:
	a1 = (force(body) ./ m)
	v1 = v0 .+ (a1 .* dt)

	# Berechne Vorhersage nach RK2:
	xnew = x0 .+ (dt/2) .* (v0 .+ v1)
	vnew = v0 .+ (dt/2) .* (a0 .+ a1)

	# Schreibe neue Werte in body:
	body[1] = xnew
	body[2] = vnew

	# Gebe Ergebnis zurück:
	return body
end

# Schnellere in situ Variante:
function fast_integ_rk2!(body, force, dt)
	dt_2 = dt * 0.5
	x0 = body[1]
	v0 = body[2]
	a0 = (force(body) ./ body[3])
	body[1] = body[1] .+ (body[2] .* dt)
	body[1] = x0 .+ dt_2 .* (v0 .+ (body[2] .+ ((force(body) ./ body[3]) .* dt)))
	body[2] = v0 .+ dt_2 .* (a0 .+ (force(body) ./ body[3]))
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
	a1  = (force(body) ./ m)							# Berechne a1 über die Kraft: f/m = a

	# Ermittel nun K3:
	x2 = x0 .+ (v1 .* dt_2)								# x2 = x0 + v1*(dt/2)
	v2 = v0 .+ (a1 .* dt_2)								# Berechne v2 via: v2 = v0 + a1*(dt/2)
	body[1] = x2										# Schreibe x2 in body, damit a2 berechnet werden kann
	a2 = (force(body) ./ m)								# Berechne a2 über die Kraft: f/m = a

	# Ermittel nun K4:
	x3 = x0 .+ (v2 .* dt)								# Berechne x3 via: x3 = x0 + v2*dt
	v3 = v0 .+ (a2 .* dt)								# Berechne v3 via: v3 = v0 + a2*dt
	body[1] = x3										# Schreibe x3 in body, damit a3 berechnet werden kann
	a3 = (force(body) ./ m)								# Berechne a3 über die Kraft: f/m = a

	# Berechne nun Vorhersage nach RK4:
	xnew = x0 .+ (dt/6) .* (v0 .+ (2 .* (v1 .+ v2)) .+ v3)
	vnew = v0 .+ (dt/6) .* (a0 .+ (2 .* (a1 .+ a2)) .+ a3)

	# Schreibe gemittelte Werte in body rein für das finale Ergebnis:
	body[1] = xnew
	body[2] = vnew

	# Gebe Ergebnis-body aus:
	return body
end

# Schnellere in situ Variante:
function fast_integ_rk4!(body, force, dt)
	dt_2 = dt * 0.5
	x0 = body[1]
	v0 = body[2]
	a0 = (force(body) ./ body[3])
	x1 = x0 .+ (v0 .* dt_2)
	v1 = v0 .+ (a0 .* dt_2)
	body[1] = x1
	a1  = (force(body) ./ body[3])
	x2 = x0 .+ (v1 .* dt_2)
	v2 = v0 .+ (a1 .* dt_2)
	body[1] = x2
	a2 = (force(body) ./ body[3])
	body[1] = x0 .+ (v2 .* dt)
	body[1] = x0 .+ (dt/6) .* (v0 .+ (2 .* (v1 .+ v2)) .+ (v0 .+ (a2 .* dt)))
	body[2] = v0 .+ (dt/6) .* (a0 .+ (2 .* (a1 .+ a2)) .+ (force(body) ./ body[3]))
end