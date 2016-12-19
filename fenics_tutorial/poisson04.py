from dolfin import *

#mesh = UnitSquareMesh(50,50)
mesh = UnitCubeMesh(10,10,10)

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

#print "assemble_system gives : (bc are applied)"
A,b = assemble_system(a,L,bc)

#print A.array()

u = Function(V)
U = u.vector() 

#"""
solver = KrylovSolver("cg","ilu")
solver.parameters["absolute_tolerance"] = 1e-7
solver.parameters["relative_tolerance"] = 1e-4
solver.parameters["maximum_iterations"] = 222
#"""
set_log_level(DEBUG)
solver.solve(A,U,b)
#solve(A,U,b)

plot(u)
plot(mesh)
interactive()


