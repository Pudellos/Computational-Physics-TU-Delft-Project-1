import numpy as np
import Lennard_Jones_potential_2D as LJ
import itertools as it

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
    def __init__(self, molecules, volume, pressure = 10**5, temperature = 273.15):
        self.molecules = molecules
        self.volume = volume # [m^3]
        self.pressure = pressure # [Pa]
        self.temperature = temperature #[K]
    
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
            potentials.append(LJ.LJ_potential(i,len(self.molecules),P=self.pressure, V=self.volume, T=self.temperature, symbol_mode=symbol_mode_gas, differenciate=differenciate_gas))
        return(potentials)
    
    def atomic_distances(self):
        '''returns distances between all iterated pairs of molecules in the simulation corresponding to pairs printed by couples_of_molecules() function'''
        distances=[]
        distances_x=[]
        distances_y=[]
        for i in range(len(self.couples_of_molecules())):
            try:
                distances.append(abs(self.couples_of_molecules()[i][0]-self.couples_of_molecules()[i][1]))
            except TypeError:
                distances_x.append(abs(self.couples_of_molecules()[i][0][0]-self.couples_of_molecules()[i][1][0]))
                distances_y.append(abs(self.couples_of_molecules()[i][0][1]-self.couples_of_molecules()[i][1][1]))
                distances.append((abs((distances_x[i])**2+(distances_y[i])**2))**(1/2))
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
            try:
                if len(chosen)==2:
                    for j in i:
                        if j[0]==chosen[0] and j[1]==chosen[1]:
                            force+=-1*temp_potentials[index]*self.atomic_distances()[index]          
            except:
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
            