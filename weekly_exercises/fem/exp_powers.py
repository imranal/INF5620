import sympy as sm

def f(x):
    return sm.exp(-x)

def rhs(f,a,b,N):
    x = sm.Symbol("x")
    v = sm.zeros(N+1,1)
    
    for i in range(0,N+1):
        v[i] = sm.integrate(f(x)*x**(2*i+1),(x,a,b))
    return v
    
def lhs(a,b,N):
    x = sm.Symbol("x")
    A = sm.zeros(N+1,N+1)
    
    for i in range(0,N+1):
        for j in range(0,N+1):
            A[i,j] = sm.integrate(x**(2*i+1)*x**(2*j+1),(x,a,b))
    return A

def plot(x,y,f,title,plt=None):
    if plt is None:
        #import matplotlib.pylab as plt
        import scitools.easyviz as plt
    plt.figure()
    plt.plot(x,y,x,f)
    plt.legend(["approximation","exact"])
    plt.title(title)
    return plt
    
    
if __name__ == "__main__":
    
    a = 0
    b = 4
    N = 4
    c = lhs(a,b,N).inv("LU")*rhs(f,a,b,N)


    import numpy as np
    x = np.linspace(a,b,101)
    y1 = c[0]*x + c[1]*x**3 + c[2]*x**5
    y2 = y1 + c[3]*x**7
    y3 = y2 + c[4]*x**9

    f = np.vectorize(f)
    f = f(x)

    plot(x,y1,f,"N=2").show()
    plot(x,y2,f,"N=3").show()
    plot(x,y3,f,"N=4").show()
    raw_input("Press enter :")

