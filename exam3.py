import numpy as np
import sys
import scitools.easyviz as plt

n = int(sys.argv[1])


T = 102.2
t = np.linspace(0,T,n+1)
dt  = t[1] - t[0]

I = 2.0
V = 3.0
w = 1
#w = w_val*(1 - (1/24.)*w_val*w_val*dt*dt)
u = np.zeros(n+1)
wdt2 = dt*w*dt*w

u[0] = I
u[1] = u[0]*(1 - 0.5*wdt2) + dt*V

for n in range(1,n):
    u[n+1] = 2*u[n] - u[n-1] - wdt2*u[n]

u_e    = np.zeros(n+2)
u_e[:] = (V/np.sqrt(w))*np.sin(np.sqrt(w)*t[:])+I*np.cos(np.sqrt(w)*t[:])

plt.plot(t,u,t,u_e,legend=["discrete","exact"])
raw_input("press enter")
