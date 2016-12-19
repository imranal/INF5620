from dolfin import *

mesh = UnitSquareMesh(3,3)

V = FunctionSpace(mesh,"CG",1)

u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]")

f = Constant(-6.0)

def boundary(x,on_boundary):
    return on_boundary
    
bc = DirichletBC(V,u0,boundary)

u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u),grad(v))*dx
L = v*f*dx

A = assemble(a)
b = assemble(L)
print "Before bcs are applied :"
print A.array()
bc.apply(A,b)
print "After bcs are applied :"
print A.array()
u = Function(V)
U = u.vector() 
solve(A,U,b)



