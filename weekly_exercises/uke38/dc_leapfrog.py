from sys import exit,argv
from scitools.std import plot,zeros,linspace
from math import sin,exp,log,cos,pi,fabs

def forward_leaper(I,a,b,T,dt):
	"""
	Solves u'(t) = -a(t)*u(t) + b(t), u(0) = I with a forward difference scheme.
	"""

	raise_type(a)
	raise_type(b)

	dt = float(dt)
	N = int(round(T/dt)) # Rounds off to the closest integer
	T = N*dt # Recalculate the T value based upon the new N value
	u = zeros(N+1)
	t = linspace(0,T,N+1)
		
	for n in range(0,N):
		u[n+1] = u[n] - dt*a(dt*n)*u[n] + dt*b(dt*n)
	return t,u

def raise_type(a):
	try: # a small test to verify that a is a function
		a_func_value = a(3.0)
	except TypeError:
		raise TypeError('a and/or b are not functions! Aborting!')
		exit(1)

if __name__ == '__main__':
	if len(argv)>1 :
		I = float(argv[1])
		a = eval(argv[2])
		b = eval(argv[3])
		T = float(argv[4])
		dt = float(argv[5])

		t,v = forward_leaper(I,a,b,T,dt)
		plot(t,v,xlabel = 't',ylabel = 'v',title = 'Leap Frog, Leap!')
		raw_input('Press any key.')
		
	
	else :
		def poly(t):
			return t
		t,v = forward_leaper(0,poly,poly,4,0.001)
		plot(t,v,xlabel = 't',ylabel = 'v',title = 'Leap Frog, Leap!')
		raw_input('Press any key.')
