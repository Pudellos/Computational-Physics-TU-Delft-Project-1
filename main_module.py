import numpy as np
import Lennard_Jones_function_handling_errors as LJ
import itertools as it

class Molecule:
    '''Class Molecule, describes a single atom or molecule with its position in space and velocity'''
    def __init__(self, position, velocity):
        self.position = position
        self.velocity = velocity
        
class Gas:
    '''Class Gas, describes collection of instances of Class Molecule'''
    def __init__(self, molecules):
        self.molecules = molecules
    
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
    
    def LJ_potentials(self, differenciate_gas=False, symbol_mode_gas=False, P_gas=10*5, V_gas=1, T_gas=273.15):
        '''Returns Lennard Johnes potentials values corresponding to couples of molecules from output of couples_of_molecules function
        
        use help(LJ.LJ_potential) for details about parameters'''
        
        potentials=[]
        for i in self.couples_of_molecules():
            potentials.append(LJ.LJ_potential(i,len(self.molecules),P=P_gas, V=V_gas, T=T_gas, symbol_mode=symbol_mode_gas, differenciate=differenciate_gas))
        return(potentials)
    
    def force_due_to_LJ(self,chosen,P=10*5,V=1,T=273.15):
        '''Returns forces on a particle due to Lennard Johnes potential corresponding to couples of molecules from output of couples_of_molecules function
        fromula:  force = - d(potential)/dr * grad (distance between two molecules)
                where grad (distance between two molecules) = abs(r_1 - r_2)
        chosen = choose position of particle for which you want to calculate the force 
        
        use help(LJ_potentials) for info on paramteres
        '''
        force=0
        index=0
        for i in self.couples_of_molecules():
            temp_potentials=self.LJ_potentials(differenciate_gas=True, symbol_mode_gas=False, P_gas=P, V_gas=V, T_gas=T)
        for i in self.couples_of_molecules():
            for j in i:
                if j==chosen:
                    force+=-1*temp_potentials[index]*abs(self.couples_of_molecules()[index][0]-self.couples_of_molecules()[index][1])
            index+=1
        return(force)

    
#EXAMPLE OF CREATING CLASS INSTANCES AND USING FUNCTIONS
# atoms = np.array((Molecule(1,2),Molecule(2,4),Molecule(3,5),Molecule(5,6)))
# g=Gas(atoms)
# print(g.positions())
# print(g.couples_of_molecules())
# print(g.LJ_potentials())
# print(g.force_due_to_LJ(1))
# print(g.force_due_to_LJ(2))