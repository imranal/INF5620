from numpy import linspace, zeros, pi, fabs, sqrt
from sys import argv, exit
from scitools.std import plot

g = -9.81 # m/s^2

def solver(v0,z0,tp,dt,mu,rho_b,A,d,C_D,V,M,R,L,p0,temp0,rho_is_zero = False):
	""" 
	This function solves the Newton's 2nd Law with respect to the velocity.
	The problem solved is for a parachute jumper. The simulation takes into account
	the following simplification(s) :
		- One Atmospheric Layer
	We solve the differential eqns involved, using Forward Euler. Since we can not
	choose an end time for the whole simulation unless the user provides the intitial
	height for the jumper, basing upon this height we estimate the time it would take
	for the jumper to fall if he/she was only affected by the gravity force. By using 
	the equation of motion, we can determine a suitable time. However, since the user 
	provides the time for the parachute to be launched, we must take this into 
	consideration as well; i.e we must make sure that the event of pulling the 
	parachute occurs well before we reach the height 0m.
		Furthermore, since the simulation may end before the estimated time, we need
	to end the simulation after height 0m has been reached. This makes physical sense.
	The program has been adjusted in such a way that a test program can attempt to 
	observe results for rho = 0. This could probably have been implemented in a much
	better way. For further tests, say terminal velocity, user should provide 
	sufficiently large height.
		To save memory usage, only the time and velocity is stored in an array before
	being returned. The other differential eqns have been solved by simply storing the
	last values only.
	"""
	
	dt    = float(dt)
	mu    = float(mu)
	rho_b = float(rho_b)
	R 	  = float(R)
 
	# Since our problem is bounded by the height, we can estimate a maximum value
	# for the time required to run the simulation :
	T = -v0/g + sqrt(4*(g**-2)*(v0**2) - 8*z0/g)/2.0
	
	# We need to make sure that the tp value provided is less than the time required
	# to reach the ground before pulling the parachute. 
	if tp > T :
		print 'Warning: The parachute will not open within the specified height'
		exit(1)

	T = T*1.5 # We adjust for drag and buoyancy (Here an arbitrary value has been chosen)
	print 'Estimated maximum time for the fall : ',T

	N  = int(round(T/dt))
	T = N*dt
	if N < 5000: # Some arbitrary value which the user should satisfy
		print 'Warning: Smaller dt value may be required.'

	v   = zeros(N+1)
	t   = linspace(0,T,N+1)

	v[0] = v0
	z    = z0
	p    = p0

	temp = temp0 - L*z0
	M_R_T = M/(R*temp*V)
	M_g_R = M*g/R
	
	rho = p*M_R_T
	if rho_is_zero:
		rho = 0 # This also implies that the we must neglect the pressure
		p   = 0 # Refer to the original equations to understand why.
	Re = rho*d*fabs(v[0])/mu

	# Here we assume that the volume and density of the body remains constant 
	# through the whole simulation. However the area of the body will change after
	# the parachute has been pulled out. As a result the drag coefficient increases.
	a_s1 = 3*pi*mu/(rho_b*V)
	b    = M_R_T/rho_b
	a_q1 = 0.5*C_D*M_R_T
	
 	pulled = False # The parachute has not been pulled yet!

	n = 0
	for n in range(0,N): # Doing ALL the subsequent operations in the for loop is bad!
		
		if z < 0 : # The simulation will end after the jumper reaches height 0m 
			print 'Ground level reached!'
			break
		if Re >= 1:  # Use the quadratic model
			v[n+1] = v[n]*(1 - dt*fabs(v[n])*p*A*a_q1) + dt*g*(p*b - 1)
		else: 		 # Use the Stokes model
			v[n+1] = v[n]*(1 - dt*a_s1*d) + dt*g*(p*b - 1)
		# Update all the relevant values for next iteration :
		z = z + dt*v[n]
		temp = temp + L*z
		p = p*(1 - dt*M_g_R/temp)
		rho = p*M_R_T
		if rho_is_zero:
			rho = 0
			p   = 0
		Re = rho*d*fabs(v[n+1])/mu
			
		if t[n] >= tp and not pulled : # Time to pull the parachute
			C_D = 2*C_D
			d   = 15*d
			A   = 20*A
			pulled = True
	return t,v

def pre_defined_values():
	"""
	Pre defined values which can be used to run the solver.
	"""

	# These are used to define the other variables :
	m = 80.0     # Mass of the jumper and parachute in kg
	R = 0.5 	 # Cross section radius of the person falling in m

	# These are values are returned :
	v0 = 0       # start velocity
	z0 = 5000    # Initial fall height
	tp = 28.0    # The person is to pull the parachute after falling tp seconds	
	dt = 0.001   # Make sure to choose a sufficiently small time step
	mu = 1.8*10**-5	# dynamic viscosity
	rho_b = 1003 # kgm^3;
	A = pi*R**2  # Cross section area of the person falling in m
	d = R/2      # Diameter of the person falling in m
	C_D = 1.2    # Coefficient of drag.
	V = m/rho_b  # Volume of the falling object
	M  = 0.029   # kg/mol
	R  = 8.314   # Nm/(mol*K)
	L  = 0.0065  # K/m
	p0 = 1*10**5 # Pa
	temp0 = 288  # K

	return v0,z0,tp,dt,mu,rho_b,A,d,C_D,V,M,R,L,p0,temp0


if __name__ == '__main__':
	v0,z0,tp,dt,mu,rho_b,A,d,C_D,V,M,R,L,p0,temp0 = pre_defined_values()
	t,v = solver(v0,z0,tp,dt,mu,rho_b,A,d,C_D,V,M,R,L,p0,temp0)
	plot(t,v,xlabel = 't',ylabel = 'v',title = 'Parachute Jumper Extreme')
	raw_input('Press any key.')
