import vertical_motion as vm
import nose.tools as nt

def test_solver():
	""" 
	In this test we neglect air resistance and buoyancy force.
	To accomplish this we can effectively set the density of air 
	to be zero. As a result we can compare the result from the solver
	with the known physical formula v = v0 - g*t.
	"""
	# This is annoying
	v0 = 0 # starts at rest
	rho = 0 # Neglect the buoyancy force as well as air resistance
	T = 8 # Choose a sufficiently large simulation
	dt = 0.1 # Apparently for smaller time steps the error increases slightly 
	g = -9.81 # Constant of gravity
	a_s = 0 # Constant to be used for the Stokes model
	a_q = 0 # Constant to be used for the quadratic model
	b = g  # Common constant for both Stokes and quadratic
	rdm = 0  # Constant used to assess the Reynolds number

	t,v = vm.solver(T,dt,a_s,a_q,b,v0,rdm)
	v_est = v0 + g*T # v = v0 + g*t
	nt.assert_almost_equal(v[-1],v_est,delta=1E-15)
	
