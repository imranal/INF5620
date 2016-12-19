class Problem:
    def __init__(self,N0=3000000,I=25000,b=0.010,d=0.0074):
        """
        The original problem is :
        
        N' = (b-d)N + I,  N(0) = N0 
        
        Introducing a fractional measure of population (for simplification) :
        
        u = N/N0
        u' = ru + f,  u(0) = 1.0
        
        Note : immigration, f, only kicks in when birth drops at some time t_r.
               See class solve. 
        """
        self.N0 = float(N0) # Initial population of the country
        self.I  = I  # Individuals immigrating to the country
        self.b  = b  # Birth rate
        self.d  = d  # Death rate
        
        # For simplification introduce fractional measure of population : u = N/N0
        self.u0 = 1.0
        self.r  = b-d
        self.f  = I/N0
        # May want to scale plot back to normal scales
        
    def define_command_line_options(self, parser=None):
        if parser==None:
            import argparse
            parser = argparse.ArgumentParser()
        
        parser.add_argument('--N0','--initial_population',type=int,default=self.N0,
                            help='initial condition : N(0)', metavar='N0')
                            
        parser.add_argument('--I','--immigrating_individuals',type=int,default=self.I, 
                            help='immigration amount : I', metavar='I')
                            
        parser.add_argument('--b','--birth_rate',type=float,default=self.b, 
                            help='birth rate b in percentage', metavar='b')
                            
        parser.add_argument('--d','--death_rate',type=float,default=self.d, 
                            help='death rate d in percentage', metavar='d')
                                              
        return parser
        

    def init_from_command_line(self, args):
        self.N0, self.I, self.b, self.d = args.N0, args.I, args.b, args.d
        
        # update these values as well :
        self.r  = self.b-self.d
        self.f  = self.I/self.N0
        
        if self.b > 1 or self.d > 1:
            raise ValueError("Birth higher or death rate given is higher than 1!")
        
        
class Solver():
    def __init__(self, problem, dt=0.1, T=101, theta=1, A=0.011, t_r=52):
        self.problem = problem
        self.dt      = float(dt)
        self.T       = T
        self.theta   = theta
        self.A  = A  # causes sudden change in birth rate at time t_r
        self.t_r = t_r

    def define_command_line_options(self, parser=None):
        parser.add_argument('--dt','--mesh_size',type=float,default=self.dt,
                                    help='define mesh size', metavar='dt')
                            
        parser.add_argument('--T','--simulation_time',type=float,default=self.T, 
                            help='define run time', metavar='T')
                            
        parser.add_argument('--A','--change_in_birth_rate_factor',type=float,  
                            default=self.A, help='define factor A greater than (b-d)',
                            metavar='A')  
                            
        parser.add_argument('--t_r','--time_when_change_in birth_rate',type=float,  
                            default=self.t_r, 
                            help='specify some time lesser than simulation time',
                            metavar='t_r')                   
     
        parser.add_argument('--theta','--finite_discretization',type=float,
                            default=self.theta,help='specify scheme : 0, 0.5 or 1',
                            metavar='theta')                   
                                              
        return parser
            
    def init_from_command_line(self, args):
        self.dt, self.T, self.theta, self.t_r, self.A = args.dt, args.T, args.theta, args.t_r , args.A
        if self.t_r > self.T:
            raise ValueError("Simulation time T=%g ends before birth rate drops "
                             "at t_r=%g." %(self.T,self.t_r))
         
    def solve(self):
        import numpy as np
        u0 = self.problem.u0
        r  = self.problem.r
        f  = self.problem.f
        dt = self.dt
        T  = self.T
        theta = self.theta
        
        N = int(round(T/dt))
        self.u = np.zeros(N+1)
        self.t = np.linspace(0,T,N+1)
        
        u = self.u
        t = self.t
        
        a = r
        f0 = 0 # immigration should be set to zero for all simulation time
        changed = False
        u[0] = u0
        for n in range(0,N):
            if (t[n] > self.t_r) and not changed:
                a  = r - self.A
                f0 = f
                changed = True
            u[n+1] = (u[n]*(1 + (1-theta)*a*dt) + dt*f0)/(1 - theta*a*dt)
        return u,t
        
        
class Visualizer:
    def __init__(self,problem,solver):
        self.problem, self.solver = problem, solver
        
    def plot(self,plt=None):
        u = self.solver.u
        t = self.solver.t
        if plt is None:
            #import scitools.std as plt
            import matplotlib.pylab as plt
        plt.plot(t,u,"r")
        plt.axvline(self.solver.t_r) # mark a straight line at t = t_r
        plt.hold("on")
        theta2name = {0:"FE", 1:"BE", 0.5:"CN"}
        name = theta2name.get(self.solver.theta,"")
        legends = ["numerical %s" %name]
        plt.legend(legends)
        plt.xlabel("t")
        plt.ylabel("u")
        plt.title("theta = %g , dt = %g" %(self.solver.theta, self.solver.dt))
        return plt
        

def main():
    problem = Problem()
    solver  = Solver(problem)
    viz     = Visualizer(problem,solver)
    
    # Read input from commandline
    parser = problem.define_command_line_options()
    parser = solver.define_command_line_options(parser)
    args   = parser.parse_args()
    problem.init_from_command_line(args)
    solver.init_from_command_line(args)
    
    # Solve and plot
    solver.solve()
    plt = viz.plot()
    plt.show()

if __name__ == "__main__":
    main()
        
        
        
        
        
        
