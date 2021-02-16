import main
import numpy as np
from main import Gas
from main import Molecule
from main import init
from main import time_evolution


g=init(2,(500,500,500))

print('positions')
print(g.positions())
print()
print('evolution happens..')
print()
g=time_evolution(g,4,0.005)
print('positions')
print(g.positions())
print()
print('calling all functions')
print()
print('velocities of particles')
print(g.velocities())
print()
print('molecules arranged into couples based on their position')
print(g.couples_of_molecules())
print()
print('atomic distances between couples of molecules')
print(g.atomic_distances())
print()
print('Lennard Jones potentials between couples of molecules')
print(g.lj_potentials())
print()
print('forces experianced by each molecule due to lennard johnes potentials')
print(g.lj_forces())
print()
print('calculation of potential energy of the simulated gas')
print(g.Ep())
print()
print('calculation of kinetic energy of the simulated gas')
print(g.Ek())
print()
print('calculation of total energy of the simulated gas')
print(g.E())