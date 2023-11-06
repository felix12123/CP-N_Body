# Es folgt nun die Implementation der Integratoren





# *****************
# ***** Euler *****
# *****************

function integ_euler(body0::Vector, acc::Function, dt::Float64)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Damit body[2] nicht 2x ausgelesen werden muss:
	v0 = body[2]

	# Berechne die neuen Orte und Geschwindigkeiten:
	v1 = v0 .+ (acc(body) .* dt)				# v1 = v + a*dt = v + (f/m)*dt
	x1 = body[1] .+ (v0 .* dt)								# x1 = x + v*dt

	# Schreibe neue Werte in body:
	body[1] = x1
	body[2] = v1

	# Gebe body zurück:
	return body
end

# Schnellere in situ Variante:
function fast_integ_euler!(body::Vector, acc::Function, dt::Float64)
	body[2] = body[2] .+ (acc(body) .* dt)
	body[1] = body[1] .+ (body[2] .* dt)
end





# ************************
# ***** Euler-Cromer *****
# ************************

function integ_euler_cromer(body0::Vector, acc::Function, dt::Float64)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Damit body[2] nicht 2x ausgelesen werden muss:
	v0 = body[2]

	# Berechne die neuen Orte und Geschwindigkeiten:
	v1 = v0 .+ (acc(body) .* dt)							# v1 = v0 + a0*dt = v + (f/m)*dt
	x1 = body[1] .+ ((1/2) .* (v0 .+ v1) .* dt)				# x1 = x0 + (1/2)*(v0 + v1)dt

	# Schreibe neue Werte in body:
	body[1] = x1
	body[2] = v1

	# Gebe body zurück:
	return body
end

# Schnellere in situ Variante:
function fast_integ_euler_cromer!(body::Vector, acc::Function, dt::Float64)
	v0 = body[2]
	body[2] = body[2] .+ (acc(body) .* dt)
	body[1] = body[1] .+ ((1/2) .* (v0 .+ body[2]) .* dt)
end





# ***************************
# ***** Velocity-Verlet *****
# ***************************

function integ_velocity_verlet(body0::Vector, acc::Function, dt::Float64)
	
	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Erforderliche Größen:
	# Speichere x0, v0 und a0:
	x0 = body[1]
	v0 = body[2]
	a0 = acc(body)

	# **********************************************************************
	# Velocity-Verlet für konstante Zeitschritte:
	#	# Berechnungsvorschrift für den neuen Ort:
	#	x1 = x0 .+ (v0 .* dt) + ((1/2) * a0 * dt^2)

	#	# Schreibe neues x (= x1) in body, um neues a (= a1) zu berechnen.
	#	# Dieses benötigt man um neues v (= v1) zu berechnen:
	#	body[1] = x1
	#	a1  = acc(body)
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
	a1  	= acc(body)
	v1 		= v1_05 .+ (1/2 .* a1 .* dt)				# (1.40)

	# Schreibe neues v (= v1) in body:
	body[2] = v1

	# Gebe Ergebnis zurück:
	return body
end

# Schnellere (Kick-Drift-Kick-Variante) in situ Variante:
function fast_integ_velocity_verlet!(body::Vector, acc::Function, dt::Float64)
	v1_05 = body[2] .+ ((1/2) .* acc(body) .* dt)
	body[1] = body[1] .+ (v1_05 .* dt)	
	body[2]	= v1_05 .+ (1/2 .* acc(body) .* dt)
end






# *************************
# ***** Runge-Kutta 2 *****
# *************************

# Hier der RK2-Integrator bzw. Heun-Integrator

function integ_rk2(body0::Vector, acc::Function, dt::Float64)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Speichere x0, v0 und a0:
	x0 = body[1]
	v0 = body[2]
	a0 = acc(body)

	# -Berechne die benötigten Werte:
	x1 = x0 .+ (v0 .* dt)
	
	# Schreibe neue Werte in body, damit neue Beschleunigung berechnet werden kann:
	body[1] = x1

	# Weiter in der Berechnung:
	a1 = acc(body)
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
function fast_integ_rk2!(body::Vector, acc::Function, dt::Float64)
	dt_2 = dt * 0.5
	x0 = body[1]
	v0 = body[2]
	a0 = acc(body)
	body[1] = body[1] .+ (body[2] .* dt)
	body[1] = x0 .+ dt_2 .* (v0 .+ (body[2] .+ (acc(body) .* dt)))
	body[2] = v0 .+ dt_2 .* (a0 .+ acc(body))
end





# *************************
# ***** Runge-Kutta 4 *****
# *************************

function integ_rk4(body0::Vector, acc::Function, dt::Float64)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Häufig verwendete Elemente sollten nicht durchweg auf's neue berechnet werden müssen:
	dt_2 = dt/2

	# Rechnungen zu K1:
	x0 = body[1]										# Ziehe Anfangswert von x
	v0 = body[2]										# Ziehe Anfangswert von v
	a0 = acc(body)										# Beschleunigung am Startpunkt

	# Rechnungen zu K2
	x1 = x0 .+ (v0 .* dt_2)								# Berechne x1 via: x1 = x0 + v0*(dt/2)
	v1 = v0 .+ (a0 .* dt_2)								# Berechne v1 via: v1 = v0 + a0*(dt/2)
	body[1] = x1										# Schreibe x1 in body, damit a1 berechnet werden kann
	a1  = acc(body)										# Berechne a1 über die Kraft: f/m = a

	# Ermittel nun K3:
	x2 = x0 .+ (v1 .* dt_2)								# x2 = x0 + v1*(dt/2)
	v2 = v0 .+ (a1 .* dt_2)								# Berechne v2 via: v2 = v0 + a1*(dt/2)
	body[1] = x2										# Schreibe x2 in body, damit a2 berechnet werden kann
	a2 = acc(body)										# Berechne a2 über die Kraft: f/m = a

	# Ermittel nun K4:
	x3 = x0 .+ (v2 .* dt)								# Berechne x3 via: x3 = x0 + v2*dt
	v3 = v0 .+ (a2 .* dt)								# Berechne v3 via: v3 = v0 + a2*dt
	body[1] = x3										# Schreibe x3 in body, damit a3 berechnet werden kann
	a3 = acc(body)										# Berechne a3 über die Kraft: f/m = a

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
function fast_integ_rk4!(body::Vector, acc::Function, dt::Float64)
	dt_2 = dt * 0.5
	x0 = body[1]
	v0 = body[2]
	a0 = acc(body)
	x1 = x0 .+ (v0 .* dt_2)
	v1 = v0 .+ (a0 .* dt_2)
	body[1] = x1
	a1  = acc(body)
	x2 = x0 .+ (v1 .* dt_2)
	v2 = v0 .+ (a1 .* dt_2)
	body[1] = x2
	a2 = acc(body)
	body[1] = x0 .+ (v2 .* dt)
	body[1] = x0 .+ (dt/6) .* (v0 .+ (2 .* (v1 .+ v2)) .+ (v0 .+ (a2 .* dt)))
	body[2] = v0 .+ (dt/6) .* (a0 .+ (2 .* (a1 .+ a2)) .+ acc(body))
end





# *******************
# ***** Hermite *****
# *******************

function integ_hermite(body0::Vector, accs::Function, dt::Float64)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Häufig verwendete Elemente sollten nicht durchweg auf's neue berechnet werden müssen:
	dt2   = dt^2
	dt3   = dt^3
	dt4   = dt^4
	dt5	  = dt^5
	# a 	  = accs[1]
	# a_dot = accs[2]

	# Rechnungen zu K1:
	x0 	   = body[1]																	# Ziehe Anfangswert von x
	v0 	   = body[2]																	# Ziehe Anfangswert von v
	a0, a0_dot = accs(body)
	
	# Berechne Schritt 1: Prediction von v und x
	v1_p = v0 .+ (a0 .* dt) .+ (0.5 .* a0_dot .* dt2)									# (1.41)
	x1_p = x0 .+ (v0 .* dt) .+ (0.5 .* a0 	 .* dt2) .+ ((1/6) .* a0_dot .* dt3)		# (1.42)
	body[1] = x1_p																		# Schreibe predicted x in body, damit neues a/a_dot berechnet werden kann
	body[2] = v1_p																		# Schreibe predicted v in body, damit neues a/a_dot berechnet werden kann
	
	# Berechne Schritt 2: Prediction der von a_p und a_dot_p
	a1_p, a1_dot_p = accs(body)
	
	# Berechne Schritt 3: Berechnung von a_2_dot und a_3_dot:
	a0_2_dot = -6 .* ((a0 .- a1_p) ./ dt2) .- 2 .* (((2 .* a0_dot) .+ a1_dot_p) ./ dt) 		# (1.46)
	a0_3_dot = 12 .* ((a0 .- a1_p) ./ dt3) .+ 6 .* (((a0_dot)   .+ a1_dot_p) ./ dt2)		# (1.47)
	
	# Berechne Schritt 4: Correction von v1 und x1:
	v1_corrected = v1_p .+ ((1/6)  .* a0_2_dot .* dt3) .+ ((1/24)  .* a0_3_dot .* dt4)	# (1.48)
	x1_corrected = x1_p .+ ((1/24) .* a0_2_dot .* dt4) .+ ((1/120) .* a0_3_dot .* dt5)	# (1.49)
	
	# Schreibe neue Werte in body:
	body[1] = x1_corrected
	body[2] = v1_corrected
	
	# Gebe Ergebnis-body aus:
	return body
end

# Schnellere in situ Variante:
function fast_integ_hermite!(body::Vector, accs::Function, dt::Float64)
	v0 	   = body[2]
	a0 	   = accs[1](body)
	a0_dot = accs[2](body)
	v1_p = v0 .+ (a0 * dt) .+ (0.5 .* a0_dot .* dt2)
	x1_p = body[1] .+ (v0 * dt) .+ (0.5 .* a0 	 .* dt2) .+ ((1/6) .* a0_dot .* dt3)
	body[1] = x1_p
	body[2] = v1_p
	a0_2_dot = -6 .* ((a0 - accs[1](body)) ./ dt2) .- 2 .* (((2*a0_dot) .+ accs[2](body)) ./ dt)
	a0_3_dot = 12 .* ((a0 - accs[1](body)) ./ dt3) .+ 6 .* (((a0_dot)   .+ accs[2](body)) ./ dt2)
	body[2] = v1_p .+ ((1/6)  .* a0_2_dot .* dt3) .+ ((1/24)  .* a0_3_dot .* dt4)
	body[1] = x1_p .+ ((1/24) .* a0_2_dot .* dt4) .+ ((1/120) .* a0_3_dot .* dt5)
end




# ******************************
# ***** Iterierter Hermite *****
# ******************************

function integ_iter_hermite(body0::Vector, accs::Function, dt::Float64)

	# Fertige Deepcopy von Body an, um aliasing zu vermeiden:
	body = deepcopy(body0)

	# Häufig verwendete Elemente sollten nicht durchweg auf's neue berechnet werden müssen:
	dt2   = dt^2
	dt3   = dt^3
	x0 	   = body[1]																						# Ziehe Anfangswert von x
	v0 	   = body[2]																						# Ziehe Anfangswert von v
	a0, a0_dot = accs(body)																					# Änderung der Beschleunigung am Startpunkt

	# Berechne Schritt 1: Prediction von v und x
	v1_p = v0 .+ (a0 .* dt) .+ (0.5 .* a0_dot .* dt2)														# (1.41)
	x1_p = x0 .+ (v0 .* dt) .+ (0.5 .* a0 	 .* dt2) .+ ((1/6) .* a0_dot .* dt3)							# (1.42)
	body[1] = x1_p																							# Schreibe predicted x in body, damit neues a/a_dot berechnet werden kann
	body[2] = v1_p																							# Schreibe predicted v in body, damit neues a/a_dot berechnet werden kann

	# Berechne Schritt 2: Prediction der von a_p und a_dot_p
	a1_p, a1_dot_p = accs(body) 																					# (1.43)

	# Berechne Schritt 3: Correction von v1 und x1:
	v1_corrected = v0 .+ ((1/2) .* (a1_p 	     .+ a0) .* dt) .+ ((1/12) .* (a1_dot_p .- a0_dot) .* dt2)	# (1.50)
	x1_corrected = x0 .+ ((1/2) .* (v1_corrected .+ v0) .* dt) .+ ((1/12) .* (a1_p 	  .- a0) 	 .* dt2)	# (1.51)

	# Schreibe neue Werte in body:
	body[1] = x1_corrected
	body[2] = v1_corrected

	# Gebe Ergebnis-body aus:
	return body
end

# Schnellere in situ Variante:
function fast_integ_iter_hermite!(body::Vector, accs::Function, dt::Float64)
	dt2   = dt^2
	x0 	   = body[1]
	v0 	   = body[2]
	a0 	   = accs[1](body)
	a0_dot = accs[2](body)
	body[2] = v0 .+ (a0 * dt) .+ (0.5 .* a0_dot .* dt2)
	body[1] = x0 .+ (v0 * dt) .+ (0.5 .* a0 	 .* dt2) .+ ((1/6) .* a0_dot .* dt3)
	a1_p 	 = accs[1](body)
	a1_dot_p = accs[2](body)
	body[2] = v0 .+ ((1/2) .* (a1_p 	     + a0) .* dt) .+ ((1/12) .* (a1_dot_p .- a0_dot) .* dt2)
	body[1] = x0 .+ ((1/2) .* (v1_corrected + v0) .* dt) .+ ((1/12) .* (a1_p 	  .- a0) 	 .* dt2)
end
