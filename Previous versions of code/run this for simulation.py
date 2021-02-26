import nico_test
import numpy as np
import matplotlib.pyplot as plt
from numpy import random
from nico_test import simulate
from nico_test import energy
from nico_test import force

box_dim = 3 #dimensions
size = 100 #grid size
N=2
init_pos = np.matrix([[52, 11, 9],[50, 10, 10]]) #2 particles that are close to eachother
#random.randint(size,size=(N,box_dim)) #assign random initial particle positions  
init_vel = np.matrix([[0, 0, 0],[0, 0, 0]])
#(random.rand(N,box_dim)-0.5) #assign random initial particle velocities
num_tsteps = 10000 #number of time steps
timestep = 0.001 #time step
x, v, dx, deltax, ULJ, EV = simulate(init_pos, init_vel, num_tsteps, timestep, box_dim, size)

plt.figure(1)
plt.plot(np.arange(0,num_tsteps),dx,'ro')
plt.title('distance between 2 particles in 3D')
plt.xlabel('number of timesteps')
plt.ylabel('distance (sigma)')

plt.figure(2)
plt.plot(np.arange(0,num_tsteps),ULJ,'ro')
plt.title('Potential energy')
plt.xlabel('number of timesteps')
plt.ylabel('energy (epsilon)')

plt.figure(3)
plt.plot(np.arange(0,num_tsteps),EV,'ro')
plt.title('Kinetic Energy')
plt.xlabel('number of timesteps')
plt.ylabel('energy (epsilon)')

plt.figure(4)
plt.plot(np.arange(0,num_tsteps),EV+ULJ,'ro')
plt.title('Sum of potential and kinetic energy')
plt.xlabel('number of timesteps')
plt.ylabel('energy (epsilon)')

box_dim = 2 #dimensions
size = 100 #grid size
N=2
init_pos = np.matrix([[5, 11],[99.99, 11]])
init_vel = np.matrix([[0, 0],[0, 0]])
num_tsteps = 10000 #number of time steps
timestep = 0.001 #time step
x, v, dx, deltax, ULJ, EV = simulate(init_pos, init_vel, num_tsteps, timestep, box_dim, size)

plt.figure(5)
plt.plot(x[:,0],x[:,1],'ro')
plt.title('particle that started at the edge x=99.99 is pulled to the other side and is now near x=0')