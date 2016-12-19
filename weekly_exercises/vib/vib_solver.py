import numpy as np

def solver(u0=0.1,T=24,dt=0.1):
    N = int(round(T/float(dt)))
    u = np.zeros(N+1)
    t = np.linspace(0,T,N+1)
    
    u[0] = u0
    u[1] = u[0] - 0.5*dt*dt*w*w*u[0] # follows from u' = 0 after inserting into scheme
    for n in range(1,N):
        u[n+1] = 2*u[n] - u[n-1] - dt*dt*w*w*u[n]
    return u,t

def exact(t,u0=0.1):
    return u0*np.cos(w*t)
w = 0.1
u,t = solver(T=320)
u_e = exact(t)

import matplotlib.pyplot as plt

plt.plot(t,u)
plt.hold("on")
plt.plot(t,u_e)
plt.xlabel("t")
plt.ylabel("u")
plt.legend(["Numerical","Exact"])
plt.show()
