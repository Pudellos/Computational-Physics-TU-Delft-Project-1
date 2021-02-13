import numpy as np
import Lennard_Jones_potential_2D as LJ
import itertools as it



def time_evolution_inside(gas,atoms,dt):
    '''evolves individual molecules of the gas
    gas: instance of class Gas
    atoms: numpy array of instances of class Molecule
    dt: float - timestep of the simulation'''
    
    if gas.dimension_2D():
        forces=gas.lj_forces()
        m=1
        temp=[]
        for i in range(len(forces)):
            v_x=((atoms[i]).velocity)[0]
            v_y=((atoms[i]).velocity)[1]
            v_x=(v_x)+((1/m)*(forces[i])*dt)
            v_y=(v_y)+((1/m)*(forces[i])*dt)
            v=(v_x,v_y)
            pos_x=((atoms[i]).position)[0]+(v_x*dt)
            pos_y=((atoms[i]).position)[1]+(v_y*dt)
            pos_x=pos_x+(dt*(v_x))
            pos_y=pos_y+(dt*(v_y))
            pos=(pos_x,pos_y)
            temp.append(Molecule(v,pos))
    else:
        forces=gas.lj_forces()
        m=6.6335209e-26 #[kg]
        for i in range(len(forces)):
            ((atoms[i]).velocity)=((atoms[i]).velocity)+((1/m)*(forces[i])*dt)
            ((atoms[i]).position)=((atoms[i]).position)+(dt*(atoms[i]).velocity)
        temp=[]
        for i in atoms: 
            temp.append(Molecule(i.position,i.velocity))

    return(temp)

def time_evolution(gas,N,dt,periodic=True):
    ''' evolves the motion of the molecules of gas in time
    gas: instance of clas Gas
    N: int = number of simulation steps
    dt: float = time step of the simulation
    periodic: bool = apply periodic boundary conditions'''
    i=0
    while i!=N:
        gas=Gas(time_evolution_inside(gas,gas.molecules,dt),gas.dimensions)
        i+=1
        if periodic==True:
            gas.periodic()
    return(gas)



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
        

        
    def dimension_2D(self):
        try:
            if len(self.positions()[0])==2:
                dim2=True
        except:
                dim2=False
        return(dim2)
     

    def positions(self):
        '''Returns positions of all gas molecules'''
        positions= [m.position for m in self.molecules]
        return(positions) 
        
    def velocities(self):
        '''Returns velocities of all gas molecules'''
        velocities = [m.velocity for m in self.molecules]
        return(velocities)
    
    def couples_of_molecules(self):
        '''Iterates through all possible couples of molecules and returns a list of positions of paired molecules '''
        positions_of_atoms=[m.position for m in self.molecules]
        couples=[]
        for couple in it.combinations(positions_of_atoms,2):
            couples.append(couple)
        return(couples)
    
    def lj_potentials(self, differenciate_gas=False, symbol_mode_gas=False, P_gas=True, V_gas=True, T_gas=True):
        '''Returns Lennard Johnes potentials values corresponding to couples of molecules from output of couples_of_molecules function
        
        use help(main_module.LJ.LJ_potential) for details about parameters'''
        
        potentials=[]
        for i in self.couples_of_molecules():
            if self.dimension_2D():
                for j in i:
                    potentials.append(LJ.LJ_potential(j,len(self.molecules),P=self.pressure, V=self.volume, T=self.temperature, symbol_mode=symbol_mode_gas, differenciate=differenciate_gas))
                
            else:
                potentials.append(LJ.LJ_potential(i,len(self.molecules),P=self.pressure, V=self.volume, T=self.temperature, symbol_mode=symbol_mode_gas, differenciate=differenciate_gas))
        return(potentials)
    
    def atomic_distances(self):
        '''returns distances between all iterated pairs of molecules in the simulation corresponding to pairs printed by couples_of_molecules() function'''
        distances=[]
        distances_x=[]
        distances_y=[]

        for i in range(len(self.couples_of_molecules())):
            if self.dimension_2D():
                distances_x.append(abs(self.couples_of_molecules()[i][0][0]-self.couples_of_molecules()[i][1][0]))
                distances_y.append(abs(self.couples_of_molecules()[i][0][1]-self.couples_of_molecules()[i][1][1]))
                distances.append((abs((distances_x[i])**2+(distances_y[i])**2))**(1/2))
            else:
                distances.append(abs(self.couples_of_molecules()[i][0]-self.couples_of_molecules()[i][1]))
        return(distances)

    
    def lj_force(self,chosen,P=True,V=True,T=True):
        '''Returns forces on a particle due to Lennard Johnes potential corresponding to couples of molecules from output of couples_of_molecules function
        fromula:  force = - d(potential)/dr * grad (distance between two molecules)
                where grad (distance between two molecules) = abs(r_1 - r_2)
        chosen = choose position of particle for which you want to calculate the force 
        
        use help(class_instance_name.LJ_potentials) for info on paramteres
        '''
        force=0
        index=0
        for i in self.couples_of_molecules():
                temp_potentials=self.lj_potentials(differenciate_gas=True, symbol_mode_gas=False, P_gas=self.pressure, V_gas=self.volume, T_gas=self.temperature)
        for i in self.couples_of_molecules():
            if self.dimension_2D():
                    for j in i:
                        if j[0]==chosen[0] and j[1]==chosen[1]:
                            force+=-1*temp_potentials[index]*self.atomic_distances()[index]          
            else:
                for j in i:
                    if j==chosen:
                        force+=-1*temp_potentials[index]*self.atomic_distances()[index]
            index+=1
        return(force)
    
    def lj_forces(self):
        '''Calculates forces on all molecules and returns them as a list'''
        forces=[]
        for i in self.molecules:
            forces.append(self.lj_force(i.position))
        return(forces)
    
    def periodic(self):
        if type(self.dimensions)==float or type(self.dimensions)==int:
            positions=self.positions()
            for i in range(len(self.molecules)):
                pos=((self.molecules[i]).position)
                if pos>=self.dimensions:
                    difference=abs(pos-self.dimensions)
                    ((self.molecules[i]).position)=difference
        else:
            positions=self.positions()
            for i in range(len(self.molecules)):
                pos_x=((self.molecules[i]).position)[0]
                pos_y=((self.molecules[i]).position)[1]
                difference_x=False
                if pos_x>=self.dimensions[0]:
                    difference_x=abs(pos_x-self.dimensions[0])
                    ((self.molecules[i]).position)=(difference_x,pos_y)

                if pos_y>=self.dimensions[1]:
                    difference_y=abs(pos_y-self.dimensions[1])
                    if difference_x:
                        ((self.molecules[i]).position)=(difference_x, difference_y)
                    else:
                        ((self.molecules[i]).position)=(pos_x, difference_y)


    def Ek(self):
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