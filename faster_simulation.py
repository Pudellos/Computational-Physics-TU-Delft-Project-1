import main
import numpy as np
from main import Gas
from main import Molecule
from main import time_evolution
from main import init

g=init(2,4,(200,200))
# atoms = np.array((Molecule(1,1),Molecule(2,2),Molecule(3,3),Molecule(4,4)))
# g=Gas(atoms,200)
print('initial')
print(g.positions())
print(g.velocities())
print(g.couples_of_molecules())
print(g.lj_potentials())
print(g.atomic_distances())
print(g.lj_forces())
print(g.E())
print()
g=time_evolution(g,11,0.000005)
print('simulation happens...')
print()
print('after simulation')
print(g.positions())
print(g.velocities())
print(g.couples_of_molecules())
print(g.lj_potentials())
print(g.atomic_distances())
print(g.lj_forces())
print(g.E())