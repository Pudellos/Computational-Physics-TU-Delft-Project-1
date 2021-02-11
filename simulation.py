import main_module
import numpy as np
from main_module import Gas
from main_module import Molecule

#specify separate molecules in the simulation:
atoms = np.array((Molecule(1,2),Molecule(2,4),Molecule(3,5),Molecule(5,6)))
#combine them into a gas
g=Gas(atoms)

#examples of class methods for various quantities:
print('positions of the gas molecules are:')
print(g.positions())
print()

print('velocities of the gas molecules are:')
print(g.velocities())
print()

print('list (listed by position) of interacting pairs of molecules in the gas:')
print(g.couples_of_molecules())
print()

print('Lennard Johnes potentials between those pairs:')
print(g.LJ_potentials())
print()

print('derivative of Lennard Johnes potentials between those pairs:')
print(g.LJ_potentials(True))
print()

print('force on a chosen molecule in the gas due to lennard johnes potentials:')
print('example: force on molecule with position = 1 due to all other molecules:')
print(g.force_due_to_LJ(1))
print('example: force on molecule with position = 2 due to all other molecules:')
print(g.force_due_to_LJ(2))

print('distances between pairs of molecules in the simulation as given by couples_of_molecules():')
print(g.atomic_distances())

print()
print('see help(function) for detailed describtion of all parameters, how to change temperature, pressure etc of the simulation')
# print(help(g.force_due_to_LJ))
