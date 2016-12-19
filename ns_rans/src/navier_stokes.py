from dolfin import *

class NSproblem :
    """
    The purpose of this code was to create a general purpose NS solver,
    later to be amended with a turbsolver/problem.
    
    This code does basically nothing except displaying the meshes for
    each subproblem or define the problems for which we have created
    the meshes.
    """
    def __init__(self,mesh=None,u0=0,p0=0,bcs=0,solver=None):
        self.mesh = mesh
        self.solver = solver
        self.u0 = u0 
        self.p0 = p0
        self.bcs = bcs
        
    def _mesh(self,hold=False):
        plot(self.mesh,interactive=hold)
    
    def defFuncSpace(self,vspace,pspace):
        self.V = VectorFunctionSpace(self.mesh,"CG",vspace)
        self.Q = FunctionSpace(self.mesh,"CG",pspace)
        self.u = TrialFunction(self.V) # u^n+1
        self.p = TrialFunction(self.Q) # p^n+1
        self.v = TestFunction(self.V)
        self.q = TestFunction(self.Q)
        self.u_1 = Function(self.V) # u^n
        self.p_1 = Function(self.Q) # p^n
    
    def getFunctions(self):
        return self.u,self.v,self.p,self.q,self.u_1,self.p_1
    
    def getSpaces(self):
        return self.V,self.Q
    
    def setPhysParams(self,nu=1e-4,rho=1.23,dt=0.001,t=0,T=1,f=None):
        self.nu  = nu 
        self.rho = rho
        self.dt  = dt
        self.t   = t
        self.T   = T
        self.f   = f # external forces
    
    def getPhysParams(self):
        return self.nu,self.rho,self.dt,self.t,self.T,self.f

class NSsolver :
    def __init__(self,problem):
        self.problem = problem

class Chorin(NSsolver):
    def __init__(self):
        NSsolver.__init__(self)       

class ChannelFlow(NSproblem): # describes a cylinder in a channel bounded by two parallel
                              # plates
    def __init__(self):
        NSproblem.__init__(self)
        self.mesh = Mesh("mesh/channel/channel.xml")
    
    def setBcs(self):
        # to understand the boundaries see the corresponding .geo file for each problem
        
        self.xmax = 4
        self.xmin = 0
        self.ymin = 0
        self.ymax = 1
        self.rcyl = 0.15 # with radius at (x,y) = (0.5,0.5)
        
        def bnd_top(x,on_boundary):
            return on_boundary and near(x[1],self.ymax)
        def bnd_bottom(x,on_boundary):
            return on_boundary and near(x[1],self.ymin)
        def bnd_cylinder(x,on_boundary):
            r = sqrt((x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5))
            return on_boundary and (r<=self.rcyl)
        def bnd_left(x,on_boundary):
            return on_boundary and near(x[0],self.xmin)
        def bnd_right(x,on_boundary):
            return on_boundary and near(x[0],self.xmax)

        u_top      = Constant((0,0))
        u_bottom   = Constant((0,0))
        u_left     = Constant((1,0)) # u @ inlet
        u_cylinder = Constant((0,0))
        # Velocities at boundaries and inlet  
        V = self.V
        Q = self.Q
        bc_val_top    = DirichletBC(V,u_top,bnd_top)
        bc_val_cyl    = DirichletBC(V,u_cylinder,bnd_cylinder)
        bc_val_left   = DirichletBC(V,u_left,bnd_left)
        bc_val_bottom = DirichletBC(V,u_bottom,bnd_bottom)

        bcu = [bc_val_top, bc_val_cyl, bc_val_left, bc_val_bottom]
        
        p_out = 0
        # Pressures at outlet    
        outflow = DirichletBC(Q,p_out,bnd_right)
        
        bcp = [outflow]
        self.bcu = bcu
        self.bcp = bcp
    
    def getBcs(self):
        return self.bcu,self.bcp 
    
    def setInitialCondition(self):
        self.u0 = Constant((0,0))
        self.p0 = Constant(0)
    
    def getInitialCondition(self):
        return self.u0,self.p0

class BackwardStepFlow(NSproblem): # This problem has been attempted to be solved in
                                   # the Chorin solver
    def __init__(self):
        NSproblem.__init__(self)
        self.mesh = Mesh("mesh/backward_step/backward_step.xml")
    
    def setBcs(self):
        self.xmax = 18.5
        self.xmin = 0
        self.ymin = 0
        self.ymax = 1
        self.ystep = 0.5
        self.xstep = 2.5
        
        def bnd_top(x,on_boundary):
            return on_boundary and near(x[1],self.ymax)
        def bnd_bottom(x,on_boundary):
            return on_boundary and near(x[1],self.ymin)
        def bnd_step(x,on_boundary):
            return on_boundary and (near(x[1],self.ystep) or near(x[0],self.xstep))
        def bnd_left(x,on_boundary):
            return on_boundary and near(x[0],self.xmin)
        def bnd_right(x,on_boundary):
            return on_boundary and near(x[0],self.xmax)

        V = self.V
        Q = self.Q
        u_top    = Constant((0,0))
        u_bottom = Constant((0,0))
        u_left   = Constant((1,0)) # u_inlet
        u_step   = Constant((0,0))
        # Velocities at boundaries and inlet  
        bc_val_top    = DirichletBC(V,u_top,bnd_top)
        bc_val_step   = DirichletBC(V,u_step,bnd_step)
        bc_val_left   = DirichletBC(V,u_left,bnd_left)
        bc_val_bottom = DirichletBC(V,u_bottom,bnd_bottom)

        bcu = [bc_val_top, bc_val_step, bc_val_left, bc_val_bottom]

        p_out = 0
        # Pressures at outlet    
        outflow = DirichletBC(Q,p_out,bnd_right)
        
        bcp = [outflow]
        self.bcu = bcu
        self.bcp = bcp
    
    def getBcs(self):
        return self.bcu,self.bcp 
        
    def setInitialCondition(self):
        self.u0 = Constant((0,0))
        self.p0 = Constant(0)
    
    def getInitialCondition(self):
        return self.u0,self.p0       
        
class CouetteFlow(NSproblem):

    def __init__(self):
        NSproblem.__init__(self)
        self.mesh = Mesh("mesh/couette/couette.xml")
        
    def setBcs(self):
        self.xmax = 2
        self.xmin = 0
        self.ymin = 0
        self.ymax = 1
    
        def bnd_top(x,on_boundary):
            return on_boundary and near(x[1],self.ymax)
        def bnd_bottom(x,on_boundary):
            return on_boundary and near(x[1],self.ymin)
        def bnd_right(x,on_boundary):
            return on_boundary and near(x[0],self.xmax)


        u_top    = Constant((1,0))
        u_bottom = Constant((0,0))
        # Velocities at boundaries and inlet  
        V = self.V
        Q = self.Q
        bc_val_top    = DirichletBC(V,u_top,bnd_top)
        bc_val_bottom = DirichletBC(V,u_bottom,bnd_bottom)

        bcu = [bc_val_top, bc_val_bottom]
        
        p_out = 0
        # Pressures at outlet    
        outflow = DirichletBC(Q,p_out,bnd_right)
        bcp = [outflow]
        
        self.bcu = bcu
        self.bcp = bcp
    
    def getBcs(self):
        return self.bcu,self.bcp  

    def setInitialCondition(self):
        self.u0 = Constant((0,0))
        self.p0 = Constant(0)
    
    def getInitialCondition(self):
        return self.u0,self.p0

if __name__ == "__main__":
    problem1 = ChannelFlow()
    problem1._mesh()

    problem2 = BackwardStepFlow()
    problem2._mesh()

    problem3 = CouetteFlow()
    problem3._mesh(hold=True)

