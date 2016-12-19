import skydiving as sd
import numpy as np
import nose.tools as nt

def test_linear_function():
    b = 0.47
    c = 1.3
    d = 2.3
    
    T  = 100
    dt = 0.1
    N = int(round(T/float(dt)))
    t = np.linspace(0,T,N+1)

    b = b*np.ones(N+1)
    c = c*np.ones(N+1)
    d = d*np.ones(N+1)

    v_e = c*t + d # exact solution
    a = (b-c)/(v_e*np.fabs(v_e)) # use the differential equation to determine a 
    
    problem = sd.Problem(a=a,b=b)
    solver  = sd.Solver(problem,dt,T,v0=d[0])
    solver.solve()
    diff = max(v_e - solver.v)
    nt.assert_almost_equal(diff,0,delta=1e-14)
    
def test_convergence():
    r = []
    E = []
    b_val = 0.47
    c_val = 1.3
    d_val = 2.3
    
    T  = 100
    dt = [0.5, 0.25, 0.125, 0.0625, 0.03125]
    
    for dt_val in dt:
        N = int(round(T/float(dt_val)))
        t = np.linspace(0,T,N+1)
        
        b = b_val*np.ones(N+1)
        c = c_val*np.ones(N+1)
        d = d_val*np.ones(N+1)
        
        v_e = c*t + d # exact solution
        a = (b-c)/(v_e*np.fabs(v_e)) # used the differential equation to determine a(t)
        
        problem = sd.Problem(a=a,b=b)
        solver  = sd.Solver(problem,dt_val,T,v0=d[0])
        solver.solve()
        E.append(sum(dt_val*np.sqrt((v_e - solver.v)**2)))
        
    for i in range(len(E)-1):
        r.append(np.log(E[i]/E[i+1])/np.log(dt[i]/dt[i+1]))
        
    nt.assert_almost_equal(r[-1],2,places=1)
