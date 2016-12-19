from dolfin import *
#set_log_level(DEBUG)

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

# Solution :
u = Function(V)

problem = LinearVariationalProblem(a,L,u,bc)
solver  = LinearVariationalSolver(problem)
"""
solver.parameters["linear_solver"] = "cg"
solver.parameters["preconditioner"] = "ilu"
cg_prm = solver.parameters["krylov_solver"]
cg_prm["absolute_tolerance"] = 1e-16
cg_prm["relative_tolerance"] = 1e-16
cg_prm["maximum_iterations"] = 39
"""
solver.solve()

u_e = interpolate(u0,V)

import numpy

print "Max error : ",numpy.abs(u_e.vector().array() - u.vector().array()).max()

center = (0.5,0.5)
print "numerical solution at center point : %.16f" % u(center)
print "exact     solution at center point : %.16f" % u0(center)
