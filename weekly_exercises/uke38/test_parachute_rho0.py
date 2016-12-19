import vertical_motion as vm
from math import sqrt
import nose.tools as nt

def test_parachute_rho0() :
	v0,z0,tp,dt,mu,rho_b,A,d,C_D,V = vm.pre_defined_values()
	t,v = vm.solver_extreme(v0,z0,tp,dt,mu,rho_b,A,d,C_D,V,rho_is_zero = True)

	g = -9.81
	T = -v0/g + sqrt(4*(g**-2)*(v0**2) - 8*z0/g)/2
	v_est = v0 + g*T # v = v0 - g*t
	nt.assert_almost_equal(v[-1],v_est,delta=1E-15)
