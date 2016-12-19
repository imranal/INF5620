from numpy import zeros,linspace, exp
from scitools.easyviz import plot, figure, hold, legend

def cooling(T0, k, T_s, t_end, dt, theta = 0.5):
    N = int(round(t_end/float(dt))) # Total amount of data points in the mesh
    T = zeros(N+1)
    T[0] = T0 # initialize
    
    for n in range(0,N):
        T[n+1] = (T[n] - k*dt*((1-theta)*T[n] - T_s))/(1.0 + theta*dt*k)
        
    t = linspace(0,t_end,N+1)
    return T,t
    
def cooling2(T0, k, T_s, t_end, temp_after_death,dt, theta = 0.5): # runs untill t_end or T_s is reached
    N = int(round(t_end/float(dt))) # Total amount of data points in the mesh
    T = zeros(N+1)
    T[0] = T0 # initialize
    
    n = 0
    while n<N and T[n]>=temp_after_death-0.0001 :
        T[n+1] = (T[n] - k*dt*((1-theta)*T[n] - T_s))/(1.0 + theta*dt*k)
        n = n + 1

    t = linspace(0,t_end,N+1)
    return T,t
    
    
def exact_sol(T0,k,T_s,t_end,dt):
    N = int(round(t_end/float(dt)))
    t = linspace(0,t_end,N+1)
    exact = (T0-T_s)*exp(-t*k) + T_s
    return exact,t
    

if __name__ == "__main__":
    T,t = cooling(37.5, 0.1, 14.0, 60.0, 5, 0.5)
    plot(t,T,"r--o")
    T_exact,t = exact_sol(37.5, 0.1, 14.0, 60.0, 5)
    #figure()
    hold("on")
    plot(t,T_exact)
    legend(["discrete sol","exact sol"])

    print "maximum error :" ,max(abs(T-T_exact))
    raw_input("press enter")


