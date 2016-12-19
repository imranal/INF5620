from dolfin import *

import sys
if len(sys.argv) > 1:  
    dim = int(sys.argv[1])
    poly_order = int(sys.argv[2])
    n = int(sys.argv[3])
else:
    dim = 2
    poly_order = 1

def _Mesh(dim=2, n=50,user_defined=False,userMesh=None):
    if user_defined: 
        return userMesh()
    if dim == 1:
        return UnitIntervalMesh(n)
    elif dim == 2:
        return UnitSquareMesh(n,n)
    elif dim == 3:
        return UnitCubeMesh(20,20,20)
         

mesh = _Mesh(dim)
V    = FunctionSpace(mesh,"CG",poly_order)

u = TrialFunction(V)
v = TestFunction(V)


T  = 0.05
dt = 0.001
t  = 0
rho = 1.0
dtbyrho = Constant(dt/float(rho))


if dim == (2 or 3):
    f = Constant(0)
else:
    f = Expression("-rho*pow(x[0],3)/3 + rho*x[0]*x[0]/2 + 8*t*t*t*pow(x[0],7)/9 - 28*t*t*t*pow(x[0],6)/9 + 7*t*t*t*pow(x[0],5)/2 - 5*t*t*t*pow(x[0],4)/4 + 2*t*x[0] - t",t=t,rho=rho)

I = Expression("cos(pi*x[0])")
u0 = interpolate(I,V)

def a(u):
    if dim== (2 or 3):
        return Constant(1)
    else:
        return Constant(1) + u*u
        
u_ = Function(V)
while t<T+DOLFIN_EPS:
    a_ = inner(u,v)*dx +\
         dtbyrho*inner(a(u0)*nabla_grad(u),nabla_grad(v))*dx 
    L  = dtbyrho*inner(u0,v)*dx +\
         dtbyrho*inner(f,v)*dx
    solve(a_== L,u_)
    #plot(u_)
    u0.assign(u_)
    t = t + dt
    f.t = t

plot(u_,interactive=True) 
# Following was to be used for the 2D case for convergence testing  
if dim == 2:
    u_e = Expression("exp(-pi*pi*t)*cos(pi*x[0])",t=T)
    u_e = interpolate(u_e,V)
    e = u_e.vector().array() - u_.vector().array()
    import numpy

    E = numpy.sqrt(numpy.sum(e*e)/float(len(u_.vector().array())))
    print E  
    
    
    
    
    
    
