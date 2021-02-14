import main_module_2D
import numpy as np
from main_module_2D import Gas
from main_module_2D import Molecule
from main_module_2D import time_evolution
from main_module_2D import init

print('note that you can state the temperature and the pressure of the simulation. \
You can change dimensions of the box container of the simulation. You can choose to apply preiodic\
boundary conditions (automatically applied, change settings to turn it off)). You can choose any inital \
position and velocity of any number of particles (1D or 2D). You can adjust number of steps of the simulation and the timestep')

############################################################ 2D all functions called ###############################################
atoms = np.array((Molecule((1,1),(0,1)),Molecule((1,2),(1,0)),Molecule((2,1),(1,1))))
g=Gas(atoms,(2,2))
print('examples of calling all funcitons of Gas class, 2D')
print('positions of particles')
print(g.positions())
print('velocities of particles')
print(g.velocities())
print('molecules arranged into couples based on their position')
print(g.couples_of_molecules())
print('atomic distances between couples of molecules')
print(g.atomic_distances())
print('Lennard Jones potentials between couples of molecules')
print(g.lj_potentials())
print('forces experianced by each molecule due to lennard johnes potentials')
print(g.lj_forces())
print('calculation of potential energy of the simulated gas')
print(g.Ep())
print('calculation of kinetic energy of the simulated gas')
print(g.Ek())
print('calculation of total energy of the simulated gas')
print(g.E())
print()
############################################################ 1D all functions called ###############################################
atoms = np.array((Molecule(1,1),Molecule(0.5,0),Molecule(0.2,0.2)))
g=Gas(atoms,1.5)
print('examples of calling all funcitons of Gas class, 1D')
print('positions of particles')
print(g.positions())
print('velocities of particles')
print(g.velocities())
print('molecules arranged into couples based on their position')
print(g.couples_of_molecules())
print('atomic distances between couples of molecules')
print(g.atomic_distances())
print('Lennard Jones potentials between couples of molecules')
print(g.lj_potentials())
print('forces experianced by each molecule due to lennard johnes potentials')
print(g.lj_forces())
print('calculation of potential energy of the simulated gas')
print(g.Ep())
print('calculation of kinetic energy of the simulated gas')
print(g.Ek())
print('calculation of total energy of the simulated gas')
print(g.E())
############################################################ time evolution of the simulation ###############################################
print()
print("now let's evolve the system in time")
atoms = np.array((Molecule((1,1),(0,1)),Molecule((1,2),(1,0)),Molecule((1.5,1.2),(1,1))))
g=Gas(atoms,(2,2))
print('before time evolution')
print()
print('positions of particles')
print(g.positions())
print('velocities of particles')
print(g.velocities())
print('total energy')
print(g.E())
print()
E_before=g.E()

print('evolution happens')
g=time_evolution(g,50,0.05)

print()
print('after time evolution')
print()
print('positions of particles')
print(g.positions())
print('velocities of particles')
print(g.velocities())
print('total energy')
print(g.E())
E_after=g.E()

print()
print('is energy conserved?')
print('checking...')
if E_before==E_after:
    print('yes, energy is conserved!')
else:
    print('no, energy is not conserved')
    
############################################################ generate gaussian distribution of positions and velocities ###############################################
g=init(2,2,(1000,1000))
print()
print("now let's evolve the system in time")
g=init(2,2,(1000,1000))
print('before time evolution')
print()
print('positions of particles')
print(g.positions())
print('velocities of particles')
print(g.velocities())
print('total energy')
print(g.E())
print()
E_before=g.E()

print('evolution happens')
g=time_evolution(g,50,0.05,True)

print()
print('after time evolution')
print()
print('positions of particles')
print(g.positions())
print('velocities of particles')
print(g.velocities())
print('total energy')
print(g.E())
E_after=g.E()

print()
print('is energy conserved?')
print('checking...')
if E_before==E_after:
    print('yes, energy is conserved!')
else:
    print('no, energy is not conserved')