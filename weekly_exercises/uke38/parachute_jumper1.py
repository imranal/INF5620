import vertical_motion as vm
from math import pi
from scitools.std import plot

def main():
	v0 = 0 # start velocity
	T = 23 # here this time is chosen for falling 5000 m when only gravity is applied
	dt = 0.01
	C_D = 1.2  # Coeffi
	rho = 0.79   # density of air (assumed constant)
	rho_b = 1003.0 # kgm^3; 
	m = 80 # Mass of the jumper and parachute in kg
	R = 0.5 	
	A = pi*R**2    # Area of the sphere         pi*r^2
	V = m/rho_b   # Volume of the falling object
	d = 2.08    # Diameter of the sphere      2*r
	mu = 1.8*10**-5	# dynamic viscosity


	g = 9.81 # Gravitational constant
	
	a_s = 3*pi*rho*d*mu/(rho_b*V) # Constant to be used for the Stokes model
	a_q = 0.5*C_D*rho*A/(rho_b*V) # Constant to be used for the quadratic model
	b   = g*(-1 + rho/rho_b)      # A common constant for both Stokes and quad. model
	rdm = rho*d/mu                # A constant related to the Reynolds number 
	
	t,v = vm.solver(T,dt,a_s,a_q,b,v0,rdm)
	plot(t,v,xlabel = 't',ylabel = 'v',title = 'Parachute Jumper')
	raw_input('Press any key.')

if __name__ == '__main__':
	main()
