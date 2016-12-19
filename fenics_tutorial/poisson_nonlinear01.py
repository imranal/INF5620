from dolfin import * 

mesh = UnitSquareMesh(33,33)
V = FunctionSpace(mesh,"CG",4)

def DirichletBC0(x,on_boundary):
    return on_boundary and abs(x[0] < 1e-14)
    
def DirichletBC1(x,on_boundary):
    return on_boundary and abs(x[0]-1) < 1e-14

bc0 = DirichletBC(V,0,DirichletBC0)
bc1 = DirichletBC(V,1,DirichletBC1)
bcs = [bc0,bc1]   

m = 2
def q(u):
    return (1+u)**m
    
u   = TrialFunction(V)
v   = TestFunction(V)
u_k = interpolate(Constant(0),V)

f   = Constant(0) 
a   = inner(q(u_k)*nabla_grad(u),nabla_grad(v))*dx
L   = f*v*dx 

# Picard Iterations
u = Function(V)
eps = 1.0
tol = 1e-9
iterations = 0
maxiter    = 23
import numpy
while eps > tol and iterations < maxiter:
    iterations += 1
    solve(a == L, u, bcs)
    diff = u.vector().array() - u_k.vector().array()
    eps = numpy.linalg.norm(diff, ord=numpy.Inf)
    print "iter=%d : norm=%g" %(iterations,eps)
    u_k.assign(u)
