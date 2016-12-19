from numpy import linspace, zeros, pi, fabs, sqrt
from sys import argv, exit
from scitools.std import plot

g = -9.81 # m/s^2

def solver(T,dt,a_s,a_q,b,v0,rdm):
	""" 
	Adaptive solver for drag on a body falling vertically through air. The program
	adjusts the model it uses for each time step. The Reynolds number is used to 
	determine whether the Stokes model or the quadratic model should be applied. 
	
	Constant used for the Stokes model :
		a_s = 3*pi*rho*d*mu/(rho_b*V)
	Constant used for the quadratic model :
		a_q = 0.5*C_D*rho*A/(rho_b*V)
	Constant used for both models :
		b = g*(-1 + rho/rho_b)
	Constant used for the Reynold's number :
		rdm = rho*d/mu 
	It would have been better to simply allow the user to specify all the constants
	so that this short notation would not require the user to precalculate them.
	"""
	dt = float(dt)
	N = int(round(T/dt)) # Rounds off to the closest integer
	T = N*dt # Recalculate the T value based upon the new N value
	
	v = zeros(N+1)
	t = linspace(0,T,N+1) # in case we wish to plot

	v[0] = v0
	Re = rdm*fabs(v0) # Reynolds number
	for n in range(0,N):
		if Re >= 1: # Use the quadratic model
			v[n+1] = (v[n] + dt*b)/(1 + a_q*dt*abs(v[n]))
		else: 		 # Use the Stokes model if Re < 1
			v[n+1] = (v[n]*(1 - a_s*dt*0.5) + dt*b)/(1 + a_s*dt*0.5)
		Re = rdm*fabs(v[n+1]) # Update the Reynolds number

	return t,v
	
	
def read_command_line():
	if len(argv) != 11:
		print '\nSimple falling body model :'
		print 'Usage : Provide following commandline arguments'
		print '>>',argv[0],'rho rho_b d mu V A C_D v0 T dt'
		print '\trho   : Density of air [kg/m^3]'
		print '\trho_b : Density of body [kg/m^3]'
		print '\td     : Diameter of body [m]'
		print '\tmu    : Dynamic viscosity [Pa*s]'
		print '\tV     : Volume of the body [m^3]'
		print '\tA     : Area of the body [m^2]'
		print '\tC_D   : Drag coefficient'
		print '\tv0    : Initial velocity [m/s]'
		print '\tT     : Total Time [s]'
		print '\tdt    : Time step size [s]'
		
		print '\nParachute Jumper :'
		print 'Usage : Simply run the program without any commandline arguments.\n'
		
		exit(1)
	else:
		rho = float(argv[1])    # Density of air
		rho_b = float(argv[2])  # Density of body
		d =  float(argv[3])     # Diameter of body
		mu =  float(argv[4])    # Dynamic viscosity
		V =  float(argv[5])     # Volume of the body
		A = float(argv[6])      # Area of the body
		C_D = float(argv[7])    # Drag coefficient
		v0 = float(argv[8])     # Initial velocity
		T = float(argv[9])		# Total Time
		dt = float(argv[10])    # Time step size
		
		return rho,rho_b,d,mu,V,A,C_D,v0,T,dt
		

if __name__ == '__main__':
	
	if len(argv) > 1 : # Runs the simple solver
		rho,rho_b,d,mu,V,A,C_D,v0,T,dt = read_command_line()
		a_s = 3*pi*rho*d*mu/(rho_b*V) # Constant to be used for the Stokes model
		a_q = 0.5*C_D*rho*A/(rho_b*V) # Constant to be used for the quadratic model
		b = g*(-1 + rho/rho_b)        # Common constant for both Stokes and quadratic
		rdm = rho*d/mu  # Constant used to asses the Reynolds number
	
		t,v = solver(T,dt,a_s,a_q,b,v0,rdm)
		plot(t,v,xlabel = 't',ylabel = 'v',title = 'Vertical Motion')
		raw_input('Press any key.')
		
	

