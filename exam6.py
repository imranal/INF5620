from dolfin import *

mesh = UnitIntervalMesh(50)
V = FunctionSpace(mesh,"CG",2)

def left_boundary(x,on_boundary):
    return on_boundary and near(x[0],0)
    
def right_boundary(x,on_boundary):
    return on_boundary and near(x[0],1)

import sys

if len(sys.argv) < 2 :
    print "Require left and right boundary values."
    sys.exit(1)

# left bc
beta = Constant(sys.argv[1])
# right bc
gamma = Constant(sys.argv[2])

bcs = [DirichletBC(V,beta,left_boundary), DirichletBC(V,gamma,right_boundary)]

u = TrialFunction(V)
v = TestFunction(V)

a = u.dx(0)*v.dx(0)*dx
L = Constant(2)*v*dx
u = Function(V)
solve(a==L,u,bcs)

plot(u)
interactive()
