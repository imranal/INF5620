from dolfin import *
from navier_stokes import *

problem = BackwardStepFlow()
#problem = CouetteFlow()
#problem = ChannelFlow()
problem.defFuncSpace(vspace=2,pspace=1) # VxQ = P2-P1
problem.setBcs()
problem.setPhysParams()

V,Q = problem.getSpaces()
u,v,p,q,u_1,p_1 = problem.getFunctions()
bcu,bcp = problem.getBcs()
nu,rho,dt,t,T,f = problem.getPhysParams()
dt = 0.001 # default dt is set to 0.01

#problem.setInitialCondition()
#u_1,p_1 = problem.getInitialCondition()

if bcp is None:
    print "No pressure conditions have been set."

if f is None:
    f = Constant((0,0))
        
u_   = TrialFunction(V) # u* : tentative velocity
phi_ = TrialFunction(Q) # phi = p - beta*p_1
k    = Constant(dt)
beta = 0.0 # Chorin scheme
while t < T + DOLFIN_EPS:
    # Tentative velocity
    F_1 = (1/k)*(inner(u_-u_1,v))*dx + inner(grad(u_1)*u_1,v)*dx -\
          (beta/rho)*inner(p_1,div(v))*dx + nu*inner(grad(u_1),grad(v))*dx -\
          inner(f,v)*dx
    a_1 = lhs(F_1)
    L_1 = rhs(F_1)
    u_s = Function(V)
    solve(a_1 == L_1,u_s,bcu)
    # Pressure update
    # setting the normal derivative 0 
    a_2 = (dt/rho)*inner(grad(phi_),grad(q))*dx
    L_2 = -div(u_s)*q*dx
    phi_s = Function(Q)
    solve(a_2 == L_2,phi_s,bcp)
    
    du = -(dt/rho)*grad(phi_s)
    # Calculate current velocity
    a_3 = inner(u,v)*dx 
    L_3 = inner(du+u_s,v)*dx
    solve(a_3 == L_3,u_s,bcu)
    
    # Calculate current pressure
    p = (beta*p_1 + phi_s)
    
    plot(u_s)
    plot(p)
    u_1.assign(u_s)
    t += dt

