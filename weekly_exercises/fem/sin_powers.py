import exp_powers as pw
import sympy as sm

def f(x):
    return sm.sin(x)

if __name__ == "__main__":
    import numpy as np
    a = 0
    b = 5*np.pi
    N = 5

    c = pw.lhs(a,b,N).inv("LU")*pw.rhs(f,a,b,N)

    import numpy as np
    x = np.linspace(a,b,101)
    y = c[0]*x + c[1]*x**3 + c[2]*x**5 + c[3]*x**7 + c[4]*x**9 + c[5]*x**11 #+ c[6]*x**13 # +\
        #c[7]*x**15 + c[8]*x**17
    
    f = np.vectorize(f)
    f = f(x)
    
    pw.plot(x,y,f,"N=5").show()
    
    raw_input("Press enter :")

