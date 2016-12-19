class Physical_Parameters:
    def __init__(self):
        self.delta = 1e-9
        self.set_default_physical_parameters()
    
    def set_default_physical_parameters(self):
        """
        Sets default physical parameters if user does not specify them.
        """
        self.A_diver_val = 0.22
        self.A_para_val = 44
        self.C_D_diver_val = 1.2
        self.C_D_para_val = 1.8
        self.rho = 1.0
        self.rho_b = 1.0e3
        self.g = 9.81
        self.m = 100
        self.V = 0.0066
        self.t_p = 69 # pull chute
        self.t_var = 0.2*self.t_p # allows C_D and A to vary linearly in [t_p,t_p+t_var] 
                          
    def a(self,t): 
        return 0.5*self.C_D(t)*self.rho*self.A(t)/(self.rho_b*self.V)
    def b(self,t): 
        return self.g*(self.rho/self.rho_b - 1.0)
            
    def C_D(self,t):
        if t <= (self.t_p - self.delta):
            return self.C_D_diver_val
        elif (self.t_p - self.delta) < t <= (self.t_var + self.t_p - self.delta):
            c = (self.C_D_para_val-self.C_D_diver_val)/self.t_var
            d = self.C_D_diver_val
            return c*(t-self.t_p) + d
        else:
            return self.C_D_para_val
    def A(self,t):  
        if t <= (self.t_p - self.delta):
            return self.A_diver_val
        elif (self.t_p - self.delta) < t <= (self.t_var + self.t_p -self.delta):
            c = (self.A_para_val-self.A_diver_val)/self.t_var
            d = self.A_diver_val
            return c*(t-self.t_p) + d
        else:
            return self.A_para_val
               
        
    def define_command_line_options(self, parser=None):
        if parser==None:
            import argparse
            parser = argparse.ArgumentParser()

        parser.add_argument('--C_D_diver_val',default=self.C_D_diver_val,
                            metavar='Coefficient of drag; value before chute is pulled.',
                            help = 'Used in a(t)')                     
        parser.add_argument('--C_D_para_val', default=self.C_D_para_val,
                            metavar='Coefficient of drag; value after chute is pulled.', 
                            help="Used in a(t)")
        parser.add_argument('--rho',default=self.rho,metavar='density of air')
        parser.add_argument('--rho_b',default=self.rho_b,metavar='density of body')
        parser.add_argument('--A_diver_val', default=self.A_diver_val,
                            metavar='Surface area; value before chute is pulled.',
                            help='Used in a(t)')                     
        parser.add_argument('--A_para_val',default=self.A_para_val,
                            metavar='Surface; value after chute is pulled.',help="Used in a(t)")
        parser.add_argument('--V',default=self.V,metavar='Volume of body')
        parser.add_argument('--g',default=self.g,metavar='Gravitational acceleration')
        parser.add_argument('--m',default=self.m,metavar='Mass of body')
        parser.add_argument('--t_p',default=self.t_p,metavar='Time to pull the chute!')
        parser.add_argument('--t_var',default=self.t_var,metavar='Time to adjust chute')
        return parser
        

    def init_from_command_line(self, args):
        """
            Initialize physical parameters from command line.
        """
        self.rho, self.rho_b, self.m, self.g, self.V, self.t_p , self.t_var=\
        args.rho, args.rho_b, args.m, args.g, args.V, float(args.t_p), float(args.t_var)
        
        self.A_diver_val, self.A_para_val, self.C_D_diver_val, self.C_D_para_val =\
        args.A_diver_val, args.A_para_val, args.C_D_diver_val, args.C_D_para_val

class Problem:# Problem class should inherit the Physical_Parameters class. Need to fix this. 
    def __init__(self,a=0.1,b=0.05,phys_prms=None):
        """
        If user does not specify physical parameters, a and b will be assumed to be
        constant and applicable for only free fall. 
        """
        self.a = a
        self.b = b
        self.user_specified = True # if a and b have been specified by user or
                                   # problem objected without specifying phys_prms
        self.phys_prms = phys_prms
        
        if phys_prms is not None:
            # If user specifies physical constants :
            self.a = phys_prms.a
            self.b = phys_prms.b
            user_specified = False
                       
    def define_command_line_options(self, parser=None):
        """
        If user does not specify physical parameters, he has the option of simply stating 
        a and b, either as functions or constant values. 
        """
        if parser==None:
            import argparse
            parser = argparse.ArgumentParser()
        
        parser.add_argument('--a',default=self.a,metavar='numpy function, np.a(t)'
                            ' or simply a value',
                            help = "For numpy functions, specify as " 
                            "np.a(t), for user-defined a(t) function.")                     
        parser.add_argument('--b',default=self.b,metavar='numpy function, np.b(t)'
                            ' or simply a value', help="For numpy functions, specify as " 
                            "np.b(t), for user-defined b(t) function.")
        parser.add_argument('--reach_term_vel',type=float,
                            help = "Used to reach terminal velocity at value given.") 
                  
        return parser

    def init_from_command_line(self, args):
        """
        Used to intialize command line arguments for a(t) and b(t).
        """   
        self.a, self.b = args.a, args.b
        if args.reach_term_vel is not None: # over-write any a and b specified!
            self.reach_term_vel(args.reach_term_vel)
        
    def reach_term_vel(self,vel): # tries to emulate graphs from the sky diving article
        self.b = -9.81
        self.a = -self.b/(vel*vel) 
         
class Solver:
    def __init__(self, problem, dt=0.1, T=161, v0=0.0):
        self.problem = problem
        self.dt      = float(dt)
        self.T       = T
        self.v0      = v0  # initial velocity of the skydiver
        self.N = int(round(T/dt))
        import numpy as np
        self.t = np.linspace(0,T,self.N+1)
        
    def define_command_line_options(self, parser=None):
        parser.add_argument('--dt','--mesh_size',metavar='dt',
                            type=float,default=self.dt, help='constant mesh size')
                            
        parser.add_argument('--T','--simulation_time',metavar='[s]',
                            type=float,default=self.T)
                            
        parser.add_argument('--v0','--intial_velocity',type=float,metavar='[m/s]',
                            default=self.v0)                   
                                              
        return parser   
        
    def cmdline_args_to_arrays(self):
        """
        Converts user-defined arguments to arrays. Specifically a and b.
        """
        def usage(arg):
            print "Usage :"
            print "To provide a numpy function use np.f(t) for some function f."
            print "Examples : --%s 'np.sin(t)' --%s 't+1' --%s "%(arg,arg,arg)
            "'lambda t: 0.2 if t<t_p else 0.5'"
            print "Or simply provide some value val : --a val"
    
        import numpy as np
        dt,T = self.dt, self.T
        self.N = int(round(T/dt))
        N = self.N
        self.t = np.linspace(0,T,N+1)  
        t = self.t 

        try:
            # If problem.a is a constant:
            self.problem.a = float(self.problem.a)*np.ones(N+1)
        except:
            # Then problem.a must be a user-defined function 
            if type(self.problem.a) is str:
                try:
                    # if numpy function :
                    if self.problem.a.startswith("np."):
                        self.problem.a = eval(self.problem.a)
                    # if lambda function :
                    else: # Hoping user is not evil
                        self.problem.a = np.vectorize(eval(self.problem.a))
                        self.problem.a = self.problem.a(t)
                except:
                    usage("a") 
                    raise ValueError("Invalid parameter provided for --a.")
            else:
                if hasattr(self.problem.a,"__call__"):
                    self.problem.a = np.vectorize(self.problem.a)
                    self.problem.a = self.problem.a(t)   
                else:
                    usage("a")
                    raise ValueError("Invalid parameter provided for --a.")           
                  
        try: # Try same process for problem.b
            self.problem.b = float(self.problem.b)*np.ones(N+1)
        except:       
            if type(self.problem.b) is str:
                try:
                    # if numpy function :
                    if self.problem.b.startswith("np."):
                        self.problem.b = eval(self.problem.b)
                    # if lambda function :
                    else: # Hoping user is an angel
                        self.problem.b = np.vectorize(eval(self.problem.b))
                        self.problem.b = self.problem.b(t)
                except:
                    usage("b") 
                    raise ValueError("Invalid parameter provided for --b.")
            else:
                if hasattr(self.problem.b,"__call__"):
                    self.problem.b = np.vectorize(self.problem.b)
                    self.problem.b = self.problem.b(t)   
                else:
                    usage("b")
                    raise ValueError("Invalid parameter provided for --b.")
        
       
    def init_from_command_line(self, args):
        """
        Get final parameters to run solver. Also turn a(t) and b(t) to numpy arrays.
        """
        self.dt, self.T, self.v0 = args.dt, args.T, args.v0
        self.cmdline_args_to_arrays()
            
    def solve(self,z0=None):
        import numpy as np
        v0 = self.v0
        dt = self.dt
        T  = self.T
        a  = self.problem.a
        b  = self.problem.b
        N = self.N
        self.v = np.zeros(N+1)       
        v = self.v
        t = self.t
        
        v[0] = v0
        if z0 is None: # if no drop height specified
            for n in range(0,N):
                v[n+1] = (v[n] + 0.5*(b[n+1]+b[n])*dt)/(1 + 0.5*(a[n+1]+a[n])*dt*np.fabs(v[n]))
            return v,t
        else: # will also print value of velocity before parachute opens
            t_p = self.problem.phys_prms.t_p
            self.z = np.zeros(N+1)
            z = self.z
            z[0] = z0
            n = 0
            while z[n]>1e-9 and n<N:
                v[n+1] = (v[n] + 0.5*(b[n+1]+b[n])*dt)/(1 + 0.5*(a[n+1]+a[n])*dt*np.fabs(v[n]))
                z[n+1] = z[n] + dt*0.5*(v[n+1]+v[n])
                n = n+1
            if n<N:
                print "Reached ground before simulation ended."
            tp_n = int(round(t_p/dt))-1
            print "Velocity at time when the parachute was pulled was %g."%v[tp_n]
            return v,z,t
            
        
    def forces(self):
        """
            Returns None if Solver.solve() has not been called or user did not provide 
            physical parameters (but only a and b).
        """ 
        if self.problem.user_specified == False:
            print "User did not specify physical parameters. Can not calculate forces from"\
            " velocities alone. User can use Physical_Parameters class to get either"\
            " pre-defined variables or supply them through command line arguments."
            return None        
        try:
            v = self.v
        except:
            print "Solver.solve() has not been called. Can not calculate drag force."
            return None
        
        phys_prob = self.problem.phys_prms
        m = phys_prob.m
        g = phys_prob.g
        rho = phys_prob.rho
        V = phys_prob.V
        a = self.problem.a
        b = self.problem.b 
        
        import numpy as np
        N = len(self.v) - 1
        self.F_g = m*g*np.ones(N+1)
        self.F_b = rho*g*V*np.ones(N+1)
        self.F_d = -a*v*np.fabs(v)   
        
class Visualizer:
    def __init__(self,problem,solver):
        self.problem, self.solver = problem, solver
        self.show = True
        
    def velocity(self,save=None,plt=None):
        v = self.solver.v
        t = self.solver.t
        if plt is None:
            #import scitools.std as plt
            import matplotlib.pylab as plt
        plt.figure()    
        plt.plot(t,v,"r")
        plt.legend("numerical solution")
        plt.xlabel("t - time")
        plt.ylabel("v - velocity")
        plt.title("Crank-Nicolson scheme used for dt = %g" %self.solver.dt)
        if save is not None: # save to file
            plt.savefig("velocity.png")
        return plt
    def forces(self,save=None,plt=None):
        F_g = self.solver.F_g
        F_b = self.solver.F_b
        F_d = self.solver.F_d
        t   = self.solver.t
        
        if plt is None:
            #import scitools.std as plt
            import matplotlib.pylab as plt
        plt.figure()    
        plt.plot(t,F_g,"r")
        plt.hold("on")
        plt.plot(t,F_b,"b")
        plt.plot(t,F_d,"g")
        plt.hold("off")
        plt.legend(["Gravity","Buouancy","Drag"])
        plt.xlabel("t - time")
        plt.ylabel("F - force")
        plt.title("Forces calculated")
        if save is not None: # save to file
            plt.savefig("forces.png")
        return plt
    
    def define_command_line_options(self, parser=None):
        parser.add_argument('--show',metavar='show',
                            type=bool,default=self.show, help='display graphs')
        return parser 
    def init_from_command_line(self, args):
        self.show = args.show # Why does cmdline arg always set to true ?
              

def main():
    phys_prm = Physical_Parameters()
    problem = Problem(phys_prms=phys_prm)
    solver  = Solver(problem)
    viz = Visualizer(problem,solver)
    
    # Read input from commandline
    parser = phys_prm.define_command_line_options()
    parser = problem.define_command_line_options(parser)
    parser = solver.define_command_line_options(parser)
    parser = viz.define_command_line_options(parser)
    
    args   = parser.parse_args()
    
    phys_prm.init_from_command_line(args)
    problem.init_from_command_line(args)
    solver.init_from_command_line(args)
    viz.init_from_command_line(args)
    
    # Solve and plot if viz.show
    solver.solve(z0=3000)
    solver.forces()
    plt = viz.forces(save=True)
    viz.velocity(save=True,plt=plt)
    
    if viz.show:
        plt.show()

if __name__ == "__main__":
    main()
"""
This code simulates a skydiver free falling for a certain time, whereupon a parachute is
released. Quadratic drag model has been used. Geometric averaging has been applied to the
nonlinear term v|v|. Crank Nicolson scheme has been used to solve the ODE.

The ODE to be solved :
        
v'(t) = -a(t)*v*abs(v) + b(t),  v(0) = v0 # specified in Solver class.

Run examples :

>> python skydiver.py

This will use the default physical parameters and adjust a and b accordingly. The simulation
includes free fall for a certain time period, whereupon the skydiver pulls the parachute.
The default values take into consideration a linear transition between the physical values
from free fall till parachute is fully open.

>> python skydiver.py --h

This will provide user with the ability to tailor any variable : either the physical variables
or simply the problem related variables a and b. Also note that there is a possibility to
supply a and b as constant values, numpy functions, lambda functions, and regular expressions.
If user specifies a and b, then the function supplied must take into account for parachute
being pulled at time t_p. Else only free fall will occur. Here is an example for user specified
parameters for a and b :

 python skydiving.py --a "lambda t: 0.022 if t<67 else 0.26" --b "lambda t: -44 if t<67 else -16" --t_p 67 --T 155

This produces a result which corresponds well with "skydiver article". Here the parachuter
reaches a terminal velocity of about 45 m/s and after t_p = 67 releases his parachute,
whereupon the velocity decreases "very fast". Note that these a and b functions do not take
into account the "smooth" transition that usually occurs when the jumper starts at rest and
finally pulls the parachute.

>> python skydiver.py --T 155 --t_p 65 --dt 0.1

This will run the simulation over 155 seconds, pull the parachute at 65 seconds, and use 
time steps 0.1 seconds for the simulation. Plots are always plotted at the end of simulation.
To suppress these plots, user can specify argument --show False. However, this does not work
for some strange reason. Apparently supplying an argument for show always sets it to true.
"""
