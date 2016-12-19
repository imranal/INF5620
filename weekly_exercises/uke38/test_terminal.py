import vertical_motion as vm
import nose.tools as nt
from numpy import sqrt, fabs, pi

def test_solver():
	""" 
	In this test we compare the velocities from the simulation for the last 	
	time step with the terminal velocity.
	"""

	# This is incredibly tedious : (every constant....)
	v0 = -10 # start velocity (already falling)
	T = 20 # set a relatively large simulation
	dt = 0.01
	C_D = 0.45  # Sphere
	rho = 1.2   # density of air
	rho_b = 6.8 # density of body -> for V = 4.716m^3; m = 6.8*4.716 = 32.07 kg 
	A = 3.40    # Area of the sphere         pi*r^2
	V = 4.716   # Volume of the sphere   (4/3)*pi*r^3
	d = 2.08    # Diameter of the sphere      2*r
	mu = 1.8*10**-5	# dynamic viscosity
	
	g = 9.81 # Gravitational constant
	
	a_s = 3*pi*rho*d*mu/(rho_b*V) # Constant to be used for the Stokes model
	a_q = 0.5*C_D*rho*A/(rho_b*V) # Constant to be used for the quadratic model
	b = g*(-1 + rho/rho_b) # A common constant for both Stokes and quad. model
	rdm = rho*d/mu         # A constant related to the Reynolds number 
	
	t,v = vm.solver(T,dt,a_s,a_q,b,v0,rdm)

	v_terminal = -sqrt(fabs(b)/a_q)
	if v[-1]>0: # can be risky to base the terminal velocity on the simulation
		v_terminal = -1*v_terminal 
	nt.assert_almost_equal(v[-1],v_terminal,delta=1E-7)
