class GridParameters:
    """
    GridData class can be initialized with following parameters :
    
    dt : time step size, which is by default set to 0.006
    dx : grid size in x-direction, which is by default set to 0.125
    dy : grid size in y-direction, which is by default set to 0.125
    dz : grid size in z-direction, which is by default set to None (i.e assume 2D)
    Lx : Length of the x axis
    Ly : ---- || ----- y axis
    Lz : ---- || ----- z axis
    problem3D : If user wants to run the 3D Wave solver
    
    Intention of the class was to originally to contain finite difference operators,
    but this was put on hold for now.
    """
    def __init__(self,Lx=10,Ly=10,Lz=None,dt=0.1,dx=0.30,dy=0.30,dz=None,problem3D=False):
        self.dt = dt
        self.dx = dx
        self.dy = dy
        self.Lx = Lx
        self.Ly = Ly
        self.Nx = int(round(self.Lx/float(dx)))
        self.Ny = int(round(self.Ly/float(dy)))
        self.problem3D = problem3D
        
        if problem3D:
            if Lz is None:
                self.Lz = Lx
            else:
                self.Lz = Lz
            if dz is None:
                self.dz = dx
            else:
                self.dz = dz
            self.Nz = int(round(self.Lz/float(self.dz)))

class PhysicalParameters:
    """
    This class must(can) be used in order to use the Solver class.
    The purpose of the class is to simply gather all the 
    variables needed to solve the wave equation (which have not
    already been specified in the class GridParameters). The parameter 
    class permits user to recieve pre-implemented functions, either 
    in 2D or 3D, depending on gridprm provided.
    
    Parameters :
    
    b : damping coefficient used in equation 
    I : initial condition
    V : initial condition for the time derivative
    q : wave velocity squared
    f : source term
    gridprm : GridParameter object
    
    """
    def __init__(self,b=0.0,I=None,V=None,q=None,f=None,gridprm=None):
        if gridprm is None:
            # Set default gridparameters(for a 2D problem) if none provided
            gridprm = GridParameters()
            self.gridprm = gridprm
        else:
            self.gridprm = gridprm
                
        import numpy as np  
        if I is None:
            if gridprm.problem3D:
                def I(x, y, z):
                    """Gaussian peak at (Lx/2, Ly/2)."""
                    return np.exp(-0.5*(x-gridprm.Lx/2.0)**2 - 0.5*(y-gridprm.Ly/2.0)**2\
                                 - 0.5*(z-gridprm.Lz/2.0)**2)
            else:
                def I(x, y):
                    """Gaussian peak at (Lx/2, Ly/2)."""
                    return np.exp(-0.5*(x-gridprm.Lx/2.0)**2 - 0.5*(y-gridprm.Ly/2.0)**2)
        
        if V is None:
            if gridprm.problem3D:
                def V(x, y, z):
                    return 0
            else:
                def V(x, y):
                    return 0
        
        if q is None:
            if gridprm.problem3D:
                def q(x, y, z):
                    #return 1.0
                    return np.exp(-0.5*(x-gridprm.Lx/2.0)**2 - 0.5*(y-gridprm.Ly/2.0)**2)
            else:
                def q(x, y):
                    return 1.0
        if f is None:
            if gridprm.problem3D:
                def f(x, y, z, t):
                    return 0.0
            else:
                def f(x, y, t):
                    return 0.0
            
        self.I = I
        self.V = V
        self.q = q
        self.f = f
        self.b = b
            
class Solver:
    """ 
    Solver class contains the 2D wave solver for the equation :
    
    [DtDt u]_{i,j}^n + b[D2t u]_{i,j}^n 
     = 
    [Dx q Dx u]_{i,j}^n + [Dy q Dy]_{i,j}^n + f_{i,j}^n
    
    The solve method contains the implementation for various versions.
    The solver employs the following functions :
    
    NeumannBC(...), DirichletBC(...), and OutletBC()
    
    These functions set the boundary conditions at each time level.
    
    Following examples display how to the run the solver :
    
    For default values, simply run : 
    >> solver = Solver()
    >> solver.solve()
    The default values are as following :
    
    For the Solver class :
    BC="neumann",BC_val=None,physprm=None
    
    Boundary condtition set to Neumann, and the physical parameters are set to their
    respective default values. See class PhysicalParameters.
    
    For the solve method :
    
    T=10.4, version="scalar", quiet_mode=True, display=False,
    save_figs=False, filename="figure", destination=None
    
    T is the total simulation time. quiet_mode allows user to gain detailed information 
    from the solver. version can be specified to be "scalar" and "vectorized". display
    allows the user to have the plots displayed during the simulation. save_figs,
    along with filename, and destination allow the user to store the plots of the
    simulation. At end of the simulation these plots are used to create a gif file for
    animation purposes. 
    
    
    To use the version "vectorized" for a time length T = 12.2 :
    >> solver = Solver()
    >> solver.solve(T=12.2,version="vectorized")
    
    To use Dirichlet boundary conditions :
    >> solver = Solver(BC="dirichlet")
    >> solver.solve()
    
    To use open outlet boundaru condtitions :
    >> solver = Solver(BC="openoutlet")
    >> solver.solve()
    
    To use more advance features, for example user defined functions, grid parameters etc.,
    these can be provided through the physprm object :
    
    >> import numpy as np
    >> I = lambda x,y: 0 if np.fabs(y-gridprm.Ly/2.0)  > 0.001 else 1.0
    >> b = 0
    >> V = lambda x,y: 0
    >> q = lambda x,y: 1.0
    >> f = lambda x,y,t: 0
    >> gridprm = wave2D.GridParameters(Lx=1,Ly=1,dt=0.025,dx=0.025,dy=0.025)
    >> physprm = wave2D.PhysicalParameters(b=b,I=I,V=V,q=q,f=f,gridprm=gridprm)
    >> solver = wave2D.Solver(physprm=physprm)
    >> solver.solve(T=T,version=version)
    
    This script runs a 1D plug through the 2D domain.
    
    _______________________________________________________________________________
    Following tests have been conducted with cpu time results :
    
    Variables used : 
    
    Lx=10, Ly=10, dt=0.0085, dx=0.05, dy=0.05, b=0,
    I="gaussian_peak at (Lx/2,Ly/2)", V=0, q=1, f=0,
    T=14.4
    
    Version scalar     : 3084.17 seconds.
    Version vectorized :  918.91 seconds. (3xfaster than scalar)
    Version cython     :  ???
    """     
    
    def __init__(self,BC="neumann",BC_val=None,physprm=None):
        
        self.bypass_stability = False
        self.outlet = False # assume no open outlet bcs
        if BC=="neumann":
            BC = self.NeumannBC
        
        elif BC=="dirichlet":
            BC = self.DirichletBC
            if BC_val is None:
                BC_val = 0
                self.BC_val = BC_val
            else:
                self.BC_val = BC_val
        
        elif BC=="openoutlet":
            BC = self.OutletBC
            self.outlet = True
        
        self.BC = BC    
        
        if physprm is not None:
            self.physprm = physprm
        else:
            # Use default grid and physical parameters
            physprm = PhysicalParameters()
            self.physprm = physprm
    
    ###################    Begin class Visualize     #####################
    class Visualize:
        def __init__(self,destination=None,filename="figure",plt=None):
            import os
            currentdir = os.getcwd()
            if destination is None:
                destination = currentdir
            else:
                destination = currentdir + "/"+ destination
                folder = "%s/%s/"%(destination,filename)
                # clean folder
                # not implemented yet
            self.destination = destination
            self.filename = filename
            if plt is None:
                import scitools.easyviz as plt
            self.plt = plt
            self.n = 0
        def store_fig(self,tn,Ix,Iy,u_2,u_1,u):
            n = self.n
            filename = self.filename
            destination = self.destination
            plt = self.plt
            
            if n==0:
                plt.surfc(Ix,Iy,u_2[1:-1,1:-1], title='t=%g'%tn, zlim=[-1.2, 1.2],
                        colorbar=True,caxis=[-1,1],shading='flat',view=[30,35],
                        hardcopy="%s/%s/%s_%04d.png"%(destination,filename,filename,n))
            elif n==1:
                plt.surfc(Ix,Iy,u_1[1:-1,1:-1], title='t=%g'%tn, zlim=[-1.2, 1.2],
                        colorbar=True,caxis=[-1,1],shading='flat',view=[30,35],
                        hardcopy="%s/%s/%s_%04d.png"%(destination,filename,filename,n))
            else:
                plt.surfc(Ix,Iy,u[1:-1,1:-1], title='t=%g'%tn, zlim=[-1.2, 1.2],
                        colorbar=True,caxis=[-1,1],shading='flat',view=[30,35],
                        hardcopy="%s/%s/%s_%04d.png"%(destination,filename,filename,n))   
            self.n = n + 1

        def display(self,tn,Ix,Iy,u_2,u_1,u):
            n = self.n
            plt = self.plt
            
            if n==0:
                plt.surfc(Ix,Iy,u_2[1:-1,1:-1], title='t=%g'%tn, zlim=[-1.2, 1.2],
                        colorbar=True,caxis=[-1,1],shading='flat',view=[30,35])
            elif n==1:
                plt.surfc(Ix,Iy,u_1[1:-1,1:-1], title='t=%g'%tn, zlim=[-1.2, 1.2],
                        colorbar=True,caxis=[-1,1],shading='flat',view=[30,35])
            else:
                plt.surfc(Ix,Iy,u[1:-1,1:-1], title='t=%g'%tn, zlim=[-1.2, 1.2],
                        colorbar=True,caxis=[-1,1],shading='flat',view=[30,35])   
            self.n = n + 1
        
        def movie(self):
            plt = self.plt
            folder = self.destination + "/" + self.filename + "/"
            plt.movie(folder+"%s_*.png"%(self.filename)) 
    ##################  End class Visualize  ###########################
    
    def OutletBC(self,u,Ix,Iy,Iz=None,term=None):
        if term is None: # for n=0 and special case run Neumann at boundaries
            return self.NeumannBC(u,Ix,Iy,Iz)
    
        version = self.version
        if version=="scalar":
            for j in range(0,Iy[-1]+1):
                i = Ix[0] # physical grid point 0
                #u[i-1,j] = term[i,j] + u[i+1,j]
                u[i-1,j] = u[i+1,j] 
                i = Ix[-1] # physical grid point Nx+1
                u[i+1,j] = term[i,j] + u[i-1,j]
                
            for i in range(0,Ix[-1]+1):
                j = Iy[0] # physical grid point 0
                #u[i,j-1] = term[i,j] + u[i,j+1]
                u[i,j-1] = u[i,j+1]
                j = Iy[-1] # physical grid point Ny+1
                #u[i,j+1] = term[i,j] + u[i,j-1]
                u[i,j+1] = u[i,j-1]
           #u[-1,-1] = term[-1,-1] + u[-1,-3] # Some how this point is never implemented in the loops
            u[-1,-1] = u[-1,-3]
               
        elif version=="vectorized":
            i = Ix[0]
            #u[i-1,:] = term[i+1,:] + u[i+1,:]
            u[i-1,:] = u[i+1,:]
            i = Ix[-1]
            u[i+1,:] = term[i-1,:] + u[i-1,:]
            j = Iy[0]  
            #u[:,j-1] = term[:,j+1] + u[:,j+1]
            u[:,j-1] = u[:,j+1]
            j = Iy[-1]
            #u[:,j+1] = term[:,j-1] + u[:,j-1]
            u[:,j+1] = u[:,j-1]
        return u        
     
    def NeumannBC(self,u,Ix,Iy,Iz=None,term=None):
        """
        Implements du/dn = 0 along the boundaries. This is done by updating 
        the ghost values with their "reflective" counterparts (using a centered scheme).
        """
        version = self.version
        if version=="scalar":
            for j in range(0,Iy[-1]+1):
                i = Ix[0] # physical grid point 0
                u[i-1,j] = u[i+1,j] 
                i = Ix[-1] # physical grid point Nx+1
                u[i+1,j] = u[i-1,j]
                
            for i in range(0,Ix[-1]+1):
                j = Iy[0] # physical grid point 0
                u[i,j-1] = u[i,j+1]
                j = Iy[-1] # physical grid point Ny+1
                u[i,j+1] = u[i,j-1]
            u[-1,-1] = u[-1,-3] # Some how this point is never implemented in the loops
               
        elif version=="vectorized":
            if self.physprm.gridprm.problem3D :
                i = Ix[0]
                u[i-1,:,:] = u[i+1,:,:]
                i = Ix[-1]
                u[i+1,:,:] = u[i-1,:,:]
                j = Iy[0]  
                u[:,j-1,:] = u[:,j+1,:]
                j = Iy[-1]
                u[:,j+1,:] = u[:,j-1,:]
                # ? :
                """
                k = Iz[0] # physical grid point 0
                u[:,:,k-1] = u[:,:,k+1]
                k = Iz[-1] # physical grid point Nz+1
                u[:,:,k+1] = u[:,:,k-1]
                """
            else:
                i = Ix[0]
                u[i-1,:] = u[i+1,:]
                i = Ix[-1]
                u[i+1,:] = u[i-1,:]
                j = Iy[0]  
                u[:,j-1] = u[:,j+1]
                j = Iy[-1]
                u[:,j+1] = u[:,j-1]
        return u
         
    def DirichletBC(self,u,Ix,Iy,Iz=None,term=None):
        """ 
        Dirichlet boundary conditions are set based on what value BC_val is given.
        BC_val is by default set to zero in __init__ if no value is provided. 
        
        Since the solver uses ghost points, the values in the ghost cells are updated
        to the BC_val as well. Not sure if this is necessary.
        """
        BC_val = self.BC_val
        version = self.version
        
        if version=="scalar":
            for j in range(0,Iy[-1]+1):
                i = Ix[0] # physical grid point 0
                u[i,j]   = BC_val
                u[i-1,j] = BC_val # Is it necessary to do this for the ghost point ?
                i = Ix[-1] # physical grid point Nx+1
                u[i,j]   = BC_val
                u[i+1,j] = BC_val
                
            for i in range(0,Ix[-1]+1):
                j = Iy[0] # physical grid point 0
                u[i,j]   = BC_val
                u[i,j-1] = BC_val
                j = Iy[-1] # physical grid point Ny+1
                u[i,j]   = BC_val
                u[i,j+1] = BC_val
        elif version=="vectorized":
            if self.physprm.gridprm.problem3D :
                i = Ix[0]
                u[i,:,:]   = BC_val
                u[i-1,:,:] = BC_val
                i = Ix[-1]
                u[i,:,:]   = BC_val
                u[i+1,:,:] = BC_val
                j = Iy[0]  
                u[:,j,:]   = BC_val
                u[:,j-1,:] = BC_val
                j = Iy[-1]
                u[:,j,:]   = BC_val
                u[:,j+1,:] = BC_val
                """
                k = Iz[0] # physical grid point 0
                u[:,:,k-1] = BC_val
                k = Iz[-1] # physical grid point Nz+1
                u[:,:,k+1] = BC_val
                """
            else :
                i = Ix[0]
                u[i,:]   = BC_val
                u[i-1,:] = BC_val
                i = Ix[-1]
                u[i,:]   = BC_val
                u[i+1,:] = BC_val
                j = Iy[0]  
                u[:,j]   = BC_val
                u[:,j-1] = BC_val
                j = Iy[-1]
                u[:,j]   = BC_val
                u[:,j+1] = BC_val
        return u         
            
    def solve(self,T=10.4,version="scalar",quiet_mode=True,display=False,\
              save_figs=False,filename="figure",destination=None):
        """
        The solver solves the 2D wave equation with varying coefficients inside the
        spatial derivatives. By default the solver uses the scalar version of 
        the implementation. Other implementations includes a vectorized version, and
        soon to be implemented Cython version.
        
        Variables :
        
        T : End time for simulation
        
        version : version of implementation
                  Choices include : - scalar (default)
                                    - vectorized
                                    - cython
        
        quiet_mode : If True (default), solver does not display any progress message.
                     And if False, solver gives detailed information what it is doing.
        
        save_figs : If True, saves the plots at certain time steps. These figures can
                    later be preprocessed to a single gif file. If False (defualt), 
                    no figures are saved.                            
        """
        
        if save_figs or display:
            viz = self.Visualize(destination=destination,filename=filename)
            destination = viz.destination # in case None was given as destination
        
        if version=="scalar" or version=="vectorized":
            self.version = version
        else:
            import sys
            print "Only following versions have been implemented for the solver : "
            print "scalar and vectorized"
            print "versions must be case-sensitive"
            sys.exit(1)
        
    
        import numpy as np
        self.T = T
        dt = self.physprm.gridprm.dt
        Lx = self.physprm.gridprm.Lx
        Ly = self.physprm.gridprm.Ly
        dt = self.physprm.gridprm.dt
        dx = self.physprm.gridprm.dx
        dy = self.physprm.gridprm.dy
        Nx = self.physprm.gridprm.Nx
        Ny = self.physprm.gridprm.Ny
        Nt = int(round(T/float(dt)))
        t  = np.linspace(0,T,Nt+1)
        #t = np.arange(0,T+dt,dt)
        x  = np.linspace(0,Lx,Nx+3)
        y  = np.linspace(0,Ly,Ny+3)
        self.x = x
        self.y = y
        
        if self.physprm.gridprm.problem3D :
            if version=="scalar":
                print "Only vectorized version available, setting version to vectorized."
                version="vectorized"
            Lz = self.physprm.gridprm.Lz
            dz = self.physprm.gridprm.dz
            Nz = self.physprm.gridprm.Nz 
            z  = np.linspace(0,Lz,Nz+3)
            self.z = z
            
        # Set bcs (either Dirichlet BCs or Neumann BCs) :
        BC = self.BC
       
        # Define the arrays to store the three time steps :
        # time step u^{n+1} :
        
        if self.physprm.gridprm.problem3D :
            u   = np.zeros([Nx+3,Ny+3,Nz+3]) 
            self.u = u
            # time step u^{n} :
            u_1 = np.zeros([Nx+3,Ny+3,Nz+3]) 
            self.u_1 = u_1
            # time step u^{n-1} : 
            u_2 = np.zeros([Nx+3,Ny+3,Nz+3])
            self.u_2 = u_2
        else :
            u   = np.zeros([Nx+3,Ny+3]) 
            self.u = u
            # time step u^{n} :
            u_1 = np.zeros([Nx+3,Ny+3]) 
            self.u_1 = u_1
            # time step u^{n-1} : 
            u_2 = np.zeros([Nx+3,Ny+3])
            self.u_2 = u_2
        
        if version=="vectorized":
            xv = x[:,np.newaxis] # for vectorized function evaluations
            yv = y[np.newaxis,:]
            if self.physprm.gridprm.problem3D :
                import scitools.easyviz as plt
                # what do we do here ?
                xx,yy,zz = plt.ndgrid(x,y,z)
        
        # These are the real physical points on our grid :
        Ix = range(1,u.shape[0]-1)
        self.Ix = Ix
        Iy = range(1,u.shape[1]-1)
        self.Iy = Iy
        if self.physprm.gridprm.problem3D :
            Iz = range(1,u.shape[2]-1)
            self.Iz = Iz      
    
        b = self.physprm.b
        # Assuming the following are functions
        f = self.physprm.f
        I = self.physprm.I    
        V = self.physprm.V
        q = self.physprm.q
        
        if self.physprm.gridprm.problem3D :
            def q_i_half(x,y,z,x_dx):
                return 0.5*(q(x,y,z)+q(x_dx,y,z))
            def q_j_half(x,y,z,y_dy):
                return 0.5*(q(x,y,z)+q(x,y_dy,z))
            def q_k_half(x,y,z,z_dz):
                return 0.5*(q(x,y,z)+q(x,y,z_dz))
            
            q_i_phalf_array = np.zeros([Nx+3,Ny+3,Nz+3])
            q_i_mhalf_array = np.zeros([Nx+3,Ny+3,Nz+3])
            q_j_phalf_array = np.zeros([Nx+3,Ny+3,Nz+3])
            q_j_mhalf_array = np.zeros([Nx+3,Ny+3,Nz+3])
            q_k_phalf_array = np.zeros([Nx+3,Ny+3,Nz+3])
            q_k_mhalf_array = np.zeros([Nx+3,Ny+3,Nz+3])
            
            q_i_phalf_array[1:-1,1:-1,1:-1] = q_i_half(x[1:-1],y[1:-1],z[1:-1],x[2:])
            q_i_mhalf_array[1:-1,1:-1,1:-1] = q_i_half(x[1:-1],y[1:-1],z[1:-1],x[:-2])
            q_j_phalf_array[1:-1,1:-1,1:-1] = q_j_half(x[1:-1],y[1:-1],z[1:-1],y[2:])
            q_j_mhalf_array[1:-1,1:-1,1:-1] = q_j_half(x[1:-1],y[1:-1],z[1:-1],y[:-2])
            q_k_phalf_array[1:-1,1:-1,1:-1] = q_j_half(x[1:-1],y[1:-1],z[1:-1],z[2:])
            q_k_mhalf_array[1:-1,1:-1,1:-1] = q_j_half(x[1:-1],y[1:-1],z[1:-1],z[:-2])    
        
        else :  
            # The following functions will be used as well later through-out
            def q_i_half(x,y,x_dx):
                # used arithmetic mean for q(i+1/2,j) and the other corresponding terms
                return 0.5*(q(x,y)+q(x_dx,y)) # x_dx should either be x + dx or x - dx
            def q_j_half(x,y,y_dy): 
                return 0.5*(q(x,y)+q(x,y_dy)) # y_dy should either be y + dy or y - dy
            
            if self.outlet:
                q_a = np.zeros([Nx+3,Ny+3])
                C5  = np.zeros([Nx+3,Ny+3])
            q_i_phalf_array = np.zeros([Nx+3,Ny+3])
            q_i_mhalf_array = np.zeros([Nx+3,Ny+3])
            q_j_phalf_array = np.zeros([Nx+3,Ny+3])
            q_j_mhalf_array = np.zeros([Nx+3,Ny+3])
            # Pre-compute q, f and V.
            if version=="scalar":
                # Fill all q_***_arrays with their respective values. 
                if self.outlet:
                    for i in Ix:
                        for j in Iy: 
                            q_a[i,j] = q(x[i],y[j]) # used for outlet bcs
                            C5[i,j]  = (dx/dt)/np.sqrt(q_a[i,j])
                for i in Ix:
                    for j in Iy: 
                        q_i_phalf_array[i,j] = q_i_half(x[i],y[j],x[i+1])
                        q_i_mhalf_array[i,j] = q_i_half(x[i],y[j],x[i-1])
                        q_j_phalf_array[i,j] = q_j_half(x[i],y[j],y[j+1])
                        q_j_mhalf_array[i,j] = q_j_half(x[i],y[j],y[j-1])

            elif version=="vectorized": 
                q_i_phalf_array[1:-1,1:-1] = q_i_half(x[1:-1],y[1:-1],x[2:])
                q_i_mhalf_array[1:-1,1:-1] = q_i_half(x[1:-1],y[1:-1],x[:-2])
                q_j_phalf_array[1:-1,1:-1] = q_j_half(x[1:-1],y[1:-1],y[2:])
                q_j_mhalf_array[1:-1,1:-1] = q_j_half(x[1:-1],y[1:-1],y[:-2])
                if self.outlet:
                    q_a[1:-1,1:-1] = q(x[1:-1],y[1:-1]) # used for outlet bcs
                    C5[1:-1,1:-1] = (dx/dt)/np.sqrt(q_a[1:-1,1:-1])
        
        if not self.bypass_stability:        
            if self.physprm.gridprm.problem3D :
                q_val = (np.fabs(q(x[:],y[:],z[:]))).max()
                stab_coeff = (1/float(np.sqrt(q_val)))*(1/np.sqrt(1/dx**2 + 1/dy**2 + 1/dz**2))
            else:
                q_val = (np.fabs(q(x[:],y[:]))).max()
                stab_coeff = (1/float(np.sqrt(q_val)))*(1/np.sqrt(1/dx**2 + 1/dy**2))
            
            if stab_coeff < t[1] :
                import sys 
                error_msg = "Time step dt=%f exceeds the stability limit %f."
                sys.exit(error_msg%(t[1],stab_coeff))
        
        # Following constants will be used throughout the solver
        C1   = 2.0/(2.0 + b*dt)
        Cb   = b*dt*0.5 - 1.0
        C11  = dt*(1 - 0.5*b*dt)
        dtdt = dt*dt
        C2   = dtdt/(dx*dx)
        C3   = dtdt/(dy*dy)
        if self.physprm.gridprm.problem3D :
            C4 = dtdt/(dz*dz)

        ##########  Step 1 : Initialize for the zeroth time step  ##########
        if not quiet_mode:
            print "Initializing zeroth time step" 
        
        if version=="scalar":       
            for i in Ix:
                for j in Iy:
                    u_2[i,j] = I(x[i-Ix[0]],y[j-Iy[0]])
        elif version=="vectorized":
            if self.physprm.gridprm.problem3D :
                u_2[:,:] = I(xv,yv,z[:])
            else :
                u_2[:,:] = I(xv,yv)

        if save_figs : # store graphs -> create animation at end
            import time
            if self.physprm.gridprm.problem3D : 
                plt.setp(interactive=False)
                h = plt.isosurface(xx,yy,zz,u_2,-3)
                h.setp(opacity=0.5)
                plt.shading("interp")
                plt.daspect([1,1,1])
                plt.view(3)
                plt.axis("tight")
                plt.show()
                raw_input("enter")
                
            else :
                viz.store_fig(t[0],Ix,Iy,u_2,u_1,u)        
                time.sleep(1.0)
        elif display:
            import time
            viz.display(t[0],Ix,Iy,u_2,u_1,u)
            time.sleep(1.0)
        
        ############     Update ghost points for u_2     ############
        if not quiet_mode:
            print "Updating ghost points"
        
        if self.physprm.gridprm.problem3D :
            BC(self.u_2,Ix,Iy,Iz,version)
        else :
            BC(self.u_2,Ix,Iy,version)
        
        ##########  Step 2 : Implement special case for t[n]=0  ##########
        """ 
        Here we use a leap frog scheme for u_t (short notation for 
        du/dt) and insert the unknown u(x,y,-dt) into the main scheme
        to find an expression for u_1, i.e u^{n}_{i,j}. 
        """   
        if not quiet_mode: 
            print "Initializing time step dt"
                
        # Inserting u^{n=-1} = u_1 - 2*dt*V(x,y) into the main scheme and solving for u_1 : 
        # Insert n = 0 in main scheme (see while loop) and insert expression for u^-1
        if version=="scalar":
            for i in Ix:
                for j in Iy:
                    q_i_phalf = q_i_phalf_array[i,j]
                    q_i_mhalf = q_i_mhalf_array[i,j]
                    q_j_phalf = q_j_phalf_array[i,j]
                    q_j_mhalf = q_j_mhalf_array[i,j]
                    u_1[i,j] = u_2[i,j] + C11*V(x[i],y[j])\
                                        + 0.5*C2*(q_i_phalf*(u_2[i+1,j]-u_2[i,j])-\
                                              q_i_mhalf*(u_2[i,j]-u_2[i-1,j]))\
                                        + 0.5*C3*(q_j_phalf*(u_2[i,j+1]-u_2[i,j])-\
                                              q_j_mhalf*(u_2[i,j]-u_2[i,j-1]))\
                                        + 0.5*dtdt*f(x[i],y[j],t[0])
        elif version=="vectorized":
            if self.physprm.gridprm.problem3D :
                u_1[1:-1,1:-1,1:-1] = u_2[1:-1,1:-1,1:-1] + C11*V(x[1:-1],y[1:-1],z[1:-1])\
                                     + 0.5*C2*(q_i_phalf_array[1:-1,1:-1,1:-1]*\
                                               (u_2[2:,1:-1,1:-1]-u_2[1:-1,1:-1,1:-1])-\
                                               q_i_mhalf_array[1:-1,1:-1,1:-1]*\
                                               (u_2[1:-1,1:-1,1:-1]-u_2[:-2,1:-1,1:-1]))\
                                     + 0.5*C3*(q_j_phalf_array[1:-1,1:-1,1:-1]*\
                                               (u_2[1:-1,2:,1:-1]-u_2[1:-1,1:-1,1:-1])-\
                                               q_j_mhalf_array[1:-1,1:-1,1:-1]*\
                                               (u_2[1:-1,1:-1,1:-1]-u_2[1:-1,:-2,1:-1]))\
                                     + 0.5*C4*(q_k_phalf_array[1:-1,1:-1,1:-1]*\
                                               (u_2[1:-1,1:-1,2:]-u_2[1:-1,1:-1,1:-1])-\
                                               q_k_mhalf_array[1:-1,1:-1,1:-1]*\
                                               (u_2[1:-1,1:-1,1:-1]-u_2[1:-1,1:-1,:-2]))\
                                     + 0.5*dtdt*f(x[1:-1],y[1:-1],z[1:-1],t[0])
            
            else :
                u_1[1:-1,1:-1] = u_2[1:-1,1:-1] + C11*V(x[1:-1],y[1:-1])\
                               + 0.5*C2*(q_i_phalf_array[1:-1,1:-1]*\
                                        (u_2[2:,1:-1]-u_2[1:-1,1:-1])-\
                                         q_i_mhalf_array[1:-1,1:-1]*\
                                         (u_2[1:-1,1:-1]-u_2[:-2,1:-1]))\
                               + 0.5*C3*(q_j_phalf_array[1:-1,1:-1]*\
                                        (u_2[1:-1,2:]-u_2[1:-1,1:-1])-\
                                         q_j_mhalf_array[1:-1,1:-1]*\
                                        (u_2[1:-1,1:-1]-u_2[1:-1,:-2]))\
                               + 0.5*dtdt*f(x[1:-1],y[1:-1],t[0])
        
        if save_figs :
            if self.physprm.gridprm.problem3D :
                plt.setp(interactive=False)
                h = plt.isosurface(xx,yy,zz,u_2,-3)
                h.setp(opacity=0.5)
                plt.shading("interp")
                plt.daspect([1,1,1])
                plt.view(3)
                plt.axis("tight")
                plt.show()
            
            else :
                viz.store_fig(t[1],Ix,Iy,u_2,u_1,u)           
                time.sleep(0.5)
        elif display:
            viz.display(t[1],Ix,Iy,u_2,u_1,u)           
            time.sleep(0.5)
        
        ############    Update ghost points for u_1    ############
        if not quiet_mode:
            print "Updating ghost points"
        
        if self.physprm.gridprm.problem3D :
            BC(self.u_1,Ix,Iy,Iz,version)
        else :
            BC(self.u_1,Ix,Iy,version)
            
        if self.outlet:
            term = np.zeros([Nx+3,Ny+3])
     
        ############   Run scheme for all time steps   ############   
        n = 1
        if not quiet_mode:
            print "Starting scheme"
        while t[n]<self.T-1e-9:
            if not quiet_mode:
                print "Starting time iteration for t=%.3f, n=%g"%(t[n],n)
            if version=="scalar":
                for i in Ix:
                    for j in Iy: # use arithmetic mean for q(i+1/2,j) etc.
                        q_i_phalf = q_i_phalf_array[i,j]
                        q_i_mhalf = q_i_mhalf_array[i,j]
                        q_j_phalf = q_j_phalf_array[i,j]
                        q_j_mhalf = q_j_mhalf_array[i,j]
                        u[i,j] = C1*(2*u_1[i,j] + Cb*u_2[i,j]\
                                     + C2*(q_i_phalf*(u_1[i+1,j]-u_1[i,j])-\
                                           q_i_mhalf*(u_1[i,j]-u_1[i-1,j]))\
                                     + C3*(q_j_phalf*(u_1[i,j+1]-u_1[i,j])-\
                                           q_j_mhalf*(u_1[i,j]-u_1[i,j-1]))\
                                     + dtdt*f(x[i],y[j],t[n]))
            elif version=="vectorized":
                if self.physprm.gridprm.problem3D :
                    u[1:-1,1:-1,1:-1] = C1*(2*u_1[1:-1,1:-1,1:-1] \
                                           + Cb*u_2[1:-1,1:-1,1:-1]\
                                           + C2*(q_i_phalf_array[1:-1,1:-1,1:-1]*\
                                                (u_1[2:,1:-1,1:-1]-u_1[1:-1,1:-1,1:-1])-\
                                                 q_i_mhalf_array[1:-1,1:-1,1:-1]*\
                                                (u_1[1:-1,1:-1,1:-1]-u_1[:-2,1:-1,1:-1]))\
                                           + C3*(q_j_phalf_array[1:-1,1:-1,1:-1]*\
                                                (u_1[1:-1,2:,1:-1]-u_1[1:-1,1:-1,1:-1])-\
                                                 q_j_mhalf_array[1:-1,1:-1,1:-1]*\
                                                (u_1[1:-1,1:-1,1:-1]-u_1[1:-1,:-2,1:-1]))\
                                           + C4*(q_k_phalf_array[1:-1,1:-1,1:-1]*\
                                                (u_1[1:-1,1:-1,2:]-u_1[1:-1,1:-1,1:-1])-\
                                                 q_k_mhalf_array[1:-1,1:-1,1:-1]*\
                                                (u_1[1:-1,1:-1,1:-1]-u_1[1:-1,1:-1,:-2]))\
                                           + dtdt*f(x[1:-1],y[1:-1],z[1:-1],t[n]))
                
                else :
                    u[1:-1,1:-1] = C1*(2*u_1[1:-1,1:-1] + Cb*u_2[1:-1,1:-1]\
                                 + C2*(q_i_phalf_array[1:-1,1:-1]*(u_1[2:,1:-1]-u_1[1:-1,1:-1])-\
                                       q_i_mhalf_array[1:-1,1:-1]*(u_1[1:-1,1:-1]-u_1[:-2,1:-1]))\
                                 + C3*(q_j_phalf_array[1:-1,1:-1]*(u_1[1:-1,2:]-u_1[1:-1,1:-1])-\
                                       q_j_mhalf_array[1:-1,1:-1]*(u_1[1:-1,1:-1]-u_1[1:-1,:-2]))\
                                 + dtdt*f(x[1:-1],y[1:-1],t[n]))   
                                   
            ############      Update the ghost points for u     ############
            if not quiet_mode:
                print "Updating ghost points"           
            
            if self.physprm.gridprm.problem3D :
                BC(self.u,Ix,Iy,Iz,version)
            else :
                if self.outlet:
                    if version=="scalar":
                        for i in Ix:
                            for j in Iy:
                                term[i,j] = C5[i,j]*(u_2[i,j] - u[i,j])   
                    elif version=="vectorized":
                        term[1:-1,1:-1] = C5[1:-1,1:-1]*(u_2[1:-1,1:-1] - u[1:-1,1:-1])
                    BC(self.u,Ix,Iy,version,term)
                else:
                    BC(self.u,Ix,Iy,version)
                      
            ############    Update all values for next time iteration    ############
            if self.physprm.gridprm.problem3D :
                u_2[:,:,:], u_1[:,:,:] = u_1, u
            else :
                u_2[:,:], u_1[:,:] = u_1, u
            
            n += 1 # Update time level
            if save_figs :
                if self.physprm.gridprm.problem3D :
                    plt.setp(interactive=False)
                    h = plt.isosurface(xx,yy,zz,u_2,-3)
                    h.setp(opacity=0.5)
                    plt.shading("interp")
                    plt.daspect([1,1,1])
                    plt.view(3)
                    plt.axis("tight")
                    plt.show()
                else :
                    viz.store_fig(t[n],Ix,Iy,u_2,u_1,u)
            elif display:
                viz.display(t[n],Ix,Iy,u_2,u_1,u)
                
        if save_figs:
            viz.movie()
            
        return u
       
if __name__ == "__main__":
    import sys
    import time
    version = "scalar"
    physprm = None
    if len(sys.argv) > 1:
        version = sys.argv[1]
        if len(sys.argv) > 2:
            problem3D = bool(sys.argv[2]) # still untested
            gridprm = GridParameters(problem3D=problem3D)
            physprm = PhysicalParameters(gridprm=gridprm)
    
    t0 = time.clock()
    solver = Solver(BC="openoutlet",physprm=physprm)
    solver.solve(version=version,T=22.4,save_figs=True,destination="fig",filename="gaussian")
    t1 = time.clock()
    print "Time for solver.solve(%s) call took %.2f seconds."%(version,(t1-t0))
    
