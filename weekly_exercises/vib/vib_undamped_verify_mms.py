import sympy as sm
sm.init_printing(use_unicode=False,wrap_line=False,no_global=True)
V, t, I, w, dt = sm.symbols('V t I w dt')  # global symbols
f = None  # global variable for the source term in the ODE

def ode_lhs(u):
    """Return left-hand side of ODE: u'' + w**2*u.
    u is symbolic Python function of t."""
    du = sm.diff(u(t),t)
    return sm.diff(du(t),t) + u(t)*w**2
    
def u_e(t):
    """Return chosen linear exact solution."""
    #return V*t + I
    #return t*t + V*t + I
    return t**3 + t*t + V*t + I
def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """
    return (1.0/(dt*dt))*(u(t + dt) - u(t+dt) + u(t-dt))

def D2t(u,dt):
    """Return 2nd-order finite difference for u_t.
    u is a symbolic Python function of t.
    """
    return (1.0/(2*dt))*(u(t+dt) - u(t-dt))
    
    
def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
    R = u_e(t) - u(t)
    return sm.simplify(R)
    
def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""
    R = u_e(dt) - u(dt)
    return sm.simplify(R)
    
def solver(I,V,f,w,T,dt):
    import numpy as np
    N = int(round(T/float(dt)))
    u = np.zeros(N+1)
    t = np.linspace(0,T,N+1)
    u[0] = I
    u[1] = V #0.5*(f[0]-w*w*u[0])*dt*dt + u[0] + dt*V
    
    for n in range(1,N):
        u[n+1] = (f[n]-w*w*u[n])*dt*dt + 2*u[n] - u[n-1]
    return u,t
    
    
def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' % u
    print "Initial conditions u(0)=%s, u'(0)=%s:" % \
          (u(t).subs(t, 0), sm.diff(u(t), t).subs(t, 0))
    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    f = sm.simplify(ode_lhs(u))
    # Residual in discrete equations (should be 0)
    print 'residual step1', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)

def linear():
    main(lambda t: V*t + I)
    
def quadratic():
    main(lambda t: t**2 + V*t + I)
def cubic():
    main(lambda t: t**3 + t**2 + V*t + I)

if __name__ == '__main__':
    #linear()
    #quadratic()
    #cubic()
    
    import matplotlib.pyplot as plt
    import numpy as np
    I = 0.2
    V = -0.9
    w = 1.2
    f = lambda t : I + V*t
    f = np.vectorize(f)
    T = 25
    dt = 0.01 
    N = int(round(T/float(dt)))
    t = np.linspace(0,T,N)
    f = f(t)
    plt.plot(t,f)
    plt.figure()
    u,t = solver(I,V,f,w,T,dt)
    print u
    plt.plot(t,u)
    plt.xlabel("t")
    plt.ylabel("u")
    plt.show()
    
