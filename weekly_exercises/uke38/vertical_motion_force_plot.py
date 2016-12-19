import vertical_motion as vm
from numpy import pi, fabs, ones, zeros
from scitools.std import plot, hold

def plotter(v0,T,dt,C_D,rho,rho_b,A,V,d,mu):
	
	g = 9.81 # Gravitational constant
	
	a_s = 3*pi*rho*d*mu/(rho_b*V) # Constant to be used for the Stokes model
	a_q = 0.5*C_D*rho*A/(rho_b*V) # Constant to be used for the quadratic model
	b = g*(-1 + rho/rho_b) # A common constant for both Stokes and quad. model
	rdm = rho*d/mu         # A constant related to the Reynolds number 
	
	t,v = vm.solver(T,dt,a_s,a_q,b,v0,rdm)

	F_b   = rho*g*V*ones(len(v))
	F_g   = -g*rho_b*V*ones(len(v))  

	rdm = rho*d/float(mu) 
	Re = rdm*fabs(v0)
	
	F_d_s = -a_s*v
	F_d_q = -a_q*v*fabs(v) 
	
	# The following code attempts to create a force vector based on the appropriate 
    # Reynolds value :
	
	F_d = zeros(len(v))
	R_e = rdm*fabs(v0)
	for n in range(len(v)) :
		if Re < 1 :
			F_d[n] = F_d_s[n] 
		else :
			F_d[n] = F_d_q[n]
		# Update Re :
		R_e = rdm*fabs(v[n])
	plot(t,F_b,
		 t,F_g,
		 xlabel='t',ylabel='F',
		 legend=('Buouncy Force','Gravity Force'),
		 title='Forces acting on a sphere')
	hold('on')
	plot(t,F_d,legend='Stokes and Quad for spesific Re')
	
	raw_input('Press any key.')

if __name__ == '__main__':
	
	v0 = -10 # start velocity (already falling)
	T = 13 # set a relatively large simulation
	dt = 0.01
	C_D = 0.45  # Sphere
	rho = 1.2   # density of air
	rho_b = 6.8 # density of body -> for V = 4.716m^3; mass = 6.8*4.716 = 32.07 kg 
	A = 3.40    # Area of the sphere         pi*r^2
	V = 4.716   # Volume of the sphere   (4/3)*pi*r^3
	d = 2.08    # Diameter of the sphere      2*r
	mu = 1.8*10**-5	# dynamic viscosity
	
	plotter(v0,T,dt,C_D,rho,rho_b,A,V,d,mu)
