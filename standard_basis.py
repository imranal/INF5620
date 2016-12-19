import sympy as sp
import numpy as np

def phi_r(r,X,d):
    if isinstance(X, sp.Symbol):
        h = sp.Rational(1,d) # node spacing
        nodes = [2*i*h - 1 for i in range(d+1)]
    else:
        # assume X is numerical : use floats for nodes
        nodes = np.linspace(-1,1,d+1)
    
    return Lagrange_polynomial(X,r,nodes)
    
def Lagrange_polynomial(X,i,points):
    p = 1
    for k in range(len(points)):
        if k!=i:
            p *= (X- points[k])/(points[i] - points[k])
    return p
    
def basis(d=1):
    "Return the complete basis."
    X = sp.Symbol("X")
    phi = [phi_r(r,X,d) for r in range(d+1)]
    return phi

def element_matrix(phi, Omega_e, symbolic=True):
    n = len(phi)
    A_e = sp.zeros((n,n))
    X = sp.Symbol("X")
    if symbolic:
        h = sp.Symbol("h")
    else:
        h = Omega_e[1] - Omega_e[0]    
    
    detJ = h/2.0
    for r in range(n):
        for s in range(r,n):
            A_e[r,s] = sp.integrate(phi[r]*phi[s]*detJ,(X,-1,1))
            A_e[s,r] = A_e[r,s]
    return A_e

def element_vector(f, phi, Omega_e, symbolic=True):
    n = len(phi)
    b_e = sp.zeros((n,1))
    X = sp.Symbol("X")
    if symbolic:   
        h = sp.Symbol("h")
    else:
        h = Omega_e[1] - Omega_e[0]    
    x = X*h/2.0 + (Omega_e[1] + Omega_e[0])/2 # mapping
    f = f.subs("x",x)   # substitute mapping formula for f, replaces x with x(X)
    detJ = h/2.0
    for r in range(n):
        b_e[r] = sp.integrate(f*phi[r]*detJ,(X,-1,1))
    return b_e

def assemble(nodes,elements,phi,f,symbolic=True):
    N_n,N_e = len(nodes), len(elements)
    zeros = sp.zeros
    A = zeros((N_n,N_n))
    b = zeros((N_n,1))
    for e in range((N_e)):
        Omega_e = [nodes[elements[e][0]],nodes[elements[e][-1]]]
        
        A_e = element_matrix(phi, Omega_e, symbolic)
        b_e = element_vector(f,phi,Omega_e,symbolic)
        for r in range(len(elements[e])):
            for s in range(len(elements[e])):
                A[elements[e][r],elements[e][s]] += A_e[r,s]
            b[elements[e][r]] += b_e[r]
    return A,b
    
if __name__ == "__main__":
    import sys

    d = 1
    x = sp.Symbol("x")
    f = x*(1-x)
    if len(sys.argv) > 1:
        d = int(sys.argv[1])
    print "phi ="   
    phi = basis(d=d)
    print phi
    
    Omega_e = []
    if d==1:
        Omega_e = [-1,1]
    elif d==2:
        Omega_e = [-1,0,1]
    A_e = element_matrix(phi,Omega_e)
    print "A_e ="   
    print A_e
    b_e = element_vector(f,phi,Omega_e)
    print "b_e ="   
    print b_e
    
    if d == 1:
        nodes = [0,0.5,1]
        elements = [[0,1],[1,2]]
        phi = basis(d=1)
    elif d == 2:
        nodes = [0,0.25,0.5,0.75,1,1.25,1.5]
        elements = [[0,1,2],[1,2,3],[3,4,5]]
    A,b = assemble(nodes,elements,phi,f) 
    print "A = "   
    print A
    print "b = "
    print b
    print "c ="
    c = A.LUsolve(b)
    print c
