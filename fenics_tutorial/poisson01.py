from dolfin import *

#parameters.form_compiler.quadrature_degree = 2

mesh = UnitSquareMesh(6,4)

V = FunctionSpace(mesh,"CG",2)

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
#u_e = project(u0,V)
u_e = interpolate(u0,V)

prm =  parameters["krylov_solver"]
prm["absolute_tolerance"] = 1e-14
prm["relative_tolerance"] = 1e-14
prm["maximum_iterations"] = 45
set_log_level(PROGRESS)
# set_log_level(DEBUG)
solve(a==L,u,bc,
      solver_parameters= dict(linear_solver="cg", preconditioner="ilu"))

plot(u)
plot(mesh)
print errornorm(u0,u,"L2")
#print max(u_e.vector().array() - u.vector().array())

#info(parameters, True)
interactive()
