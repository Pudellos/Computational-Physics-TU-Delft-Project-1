import numpy as np
import Lennard_Jones_potential_2D as LJ
import itertools as it
import sympy as sym

def init(dim,N,dimensions,P=10**5,T=273.15):
    atoms=np.full(N,Molecule(0,0))
    g=Gas(atoms,dimensions,P,T)
    init_v_x=g.initialize_velocity(T)
    init_v_y=g.initialize_velocity(T)
    init_pos_x=g.initialize_position()
    init_pos_y=g.initialize_position()
    index=0
    for i in atoms:
        i.velocity=(init_v_x[0][index],init_v_y[0][index])
        i.position=(init_pos_x[0][index],init_pos_y[0][index])
        molecule=Molecule(i.position,i.velocity)
        atoms[index]=molecule
        index+=1
    a=Gas(atoms,dimensions,P,T)
    return(a)

def time_evolution_inside(gas,dt):
    '''evolves individual molecules of the gas
    gas: instance of class Gas
    atoms: numpy array of instances of class Molecule
    dt: float - timestep of the simulation'''
    v,pos=gas.velocities(),gas.positions()
    forces=np.array(gas.lj_forces())
    m=6.6335209e-26 #[kg]
    v=np.array([v[i]+(1/m)*forces[i]*dt for i in range(len(gas.molecules))])
    pos=np.array([pos[i]+v[i]*dt for i in range(len(gas.molecules))])
    molecules_evolved=np.array([Molecule(v[i],pos[i]) for i in range(len(gas.molecules))])
    return(Gas(molecules_evolved,gas.dimensions,gas.pressure,gas.temperature))

def time_evolution(gas,N,dt,periodic=True):
    ''' evolves the motion of the molecules of gas in time
    gas: instance of clas Gas
    N: int = number of simulation steps
    dt: float = time step of the simulation
    periodic: bool = apply periodic boundary conditions'''
    i=0
    while i!=N:
        gas=time_evolution_inside(gas,dt)
        i+=1
    return(gas)

def periodic_positions(gas):
            positions=gas.positions()
            index=np.unique(np.where(abs(positions)>=gas.dimensions)[0])
            while index.size>0:
                difference=abs(positions-gas.dimensions)
                positions[index]=difference[index]
                index=np.unique(np.where(positions>=gas.dimensions)[0])
            return(positions)


class Molecule:
    '''Class Molecule, describes a single atom or molecule with its position in space and velocity'''
    def __init__(self, position, velocity):
        self.position = position
        self.velocity = velocity
        
class Gas:
    '''Class Gas, describes collection of instances of Class Molecule
    default pressure and temperature are the stp values
    molecules : np.ndarray = array of instances of Molecule class
    volume : float = volume of simulated gas [m^3]
    pressure : float = pressure of simulated gas [Pa]
    temperature : float = temperature of simulated gas [K]
    '''

    
    def __init__(self, molecules, dimensions, pressure = 10**5, temperature = 273.15):
        self.molecules = molecules
        self.dimensions = dimensions # [m^3]
        self.pressure = pressure # [Pa]
        self.temperature = temperature #[K]
        if type(dimensions)==float or type(dimensions)==int:
            self.volume=dimensions
        else:
            self.volume=dimensions[0]*dimensions[1]
        
    def k_B(self):
        ''' Calculation of Boltzmann constant for the simulation
    
        parameters:
        P = pressure [Pa]
        V= volume [m^3]
        T = temperature [K]
        N = number of gas molecules in the simulation '''
        return((self.pressure*self.volume)/(self.temperature*len(self.molecules)))
    
    def lj_potentials_between_pairs(self, differenciate=True):

        def equation(sigma,eps,r_abs):
            return((4*eps*((sigma/r_abs)**12-(sigma/r_abs)**6)))
    
        sigma = sym.Symbol('sigma')
        eps = sym.Symbol('eps')
        r_abs=sym.Symbol('r_abs')
    
        sigma_value=3.405e-10
        eps_value=119.8/self.k_B()
        r=np.zeros(len(self.couples_of_molecules()))
        if self.dimension_2D():
            r=np.array([np.sqrt((abs(i[0][0]-i[1][0]))**2+(abs(i[0][1]-i[1][1]))**2) for i in self.couples_of_molecules()])
        potentials=np.zeros(len(self.couples_of_molecules()))
        eqn=equation(sigma,eps,r_abs)
        potentials=np.array([eqn.evalf(subs={sigma: sigma_value, eps: eps_value, r_abs:i}) for i in r])
        
        if differenciate:
            eqn=sym.diff(equation(sigma,eps,r_abs), r_abs)
            potentials=np.array([eqn.evalf(subs={sigma: sigma_value, eps: eps_value, r_abs:i}) for i in r])
        return(potentials)
    
    def lj_potentials(self):
        k=self.lj_potentials_between_pairs()
        potentials=np.zeros(len(self.molecules))
        for i in range(len(self.molecules)):
            chosen=self.positions()[i]
            if chosen in self.couples_of_molecules():
                index=np.unique(np.where(self.couples_of_molecules()==chosen)[0])
                potentials[i]=sum(k[index])
        return(potentials)

    def dimension_2D(self):
        try:
            if len(self.positions()[0])==2:
                dim2=True
        except:
                dim2=False
        return(dim2)
     

    def positions(self):
        '''Returns positions of all gas molecules'''
        positions= np.array([m.position for m in self.molecules])
        return(positions) 
        
    def position_periodic(self):
        return(periodic_positions(self))
        
    def velocities(self):
        '''Returns velocities of all gas molecules'''
        velocities = np.array([m.velocity for m in self.molecules])
        return(velocities)
    
    def couples_of_molecules(self):
        '''Iterates through all possible couples of molecules and returns a list of positions of paired molecules '''
        couples=np.zeros(len(self.positions()))
        couples=np.array([i for i in it.combinations(self.positions(),2)])
        return(couples)
    
    def atomic_distances(self):
        '''returns distances between all iterated pairs of molecules in the simulation corresponding to pairs printed by couples_of_molecules() function'''
        r=np.zeros(len(self.couples_of_molecules()))
        if self.dimension_2D():
            r=np.array([np.sqrt((abs(i[0][0]-i[1][0]))**2+(abs(i[0][1]-i[1][1]))**2) for i in self.couples_of_molecules()])
        return(r)

    def lj_forces(self):
        forces=np.zeros(len(self.molecules))
        forces=np.array([i*sum(self.atomic_distances()) for i in self.lj_potentials()])
        return(forces)

    def Ek(self,T=True):
        '''returns kinetic energy of the gas in Joules'''
        M=39.948/1000 #[kg/mol]
        T=self.temperature
        R=8.314 #J/(mol*K)
        m=6.6335209e-26 #[kg]
        v_rms=np.sqrt(3*R*T/M)
        return((1/2)*m*v_rms**2)
    
    def Ep(self):
        sum=0
        for i in self.lj_potentials():
            sum+=i
        return(sum)
    
    def E(self):
        return(self.Ek()+self.Ep())
    
    def initialize_velocity(self,T):
        E=self.Ek(T)
        M=39.948/1000 #[kg/mol]
        R=8.314 #J/(mol*K)
        m=6.6335209e-26 #[kg]
        v_rms=np.sqrt(3*R*T/M)
        return(np.random.normal(v_rms, 10, size=(1, len(self.molecules))))
    
    def initial_rms_v(self,T):
        E=self.Ek(T)
        M=39.948/1000 #[kg/mol]
        R=8.314 #J/(mol*K)
        m=6.6335209e-26 #[kg]
        v_rms=np.sqrt(3*R*T/M)
        return(v_rms)
    
    def initialize_position(self):
        return(abs(np.random.normal(0, 0.0000001, size=(1, len(self.molecules)))))