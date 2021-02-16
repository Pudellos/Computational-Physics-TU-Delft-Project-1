import main
import numpy as np
from main import Gas
from main import Molecule
from main import init
from main import time_evolution

g=init(2,3,(500,500))
print('positions before')
print(g.positions())
print()
print('evolution happens..')
print()
g=time_evolution(g,4,0.005)
print('positions after')
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