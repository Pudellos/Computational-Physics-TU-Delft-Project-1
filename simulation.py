import skeleton_ours
from skeleton_ours import simulate, fcc_lattice, init_velocity
#from skeleton_ours import 
#from skeleton_ours import 
import itertools
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation
import time

dim = 3
num_tsteps = 20000
num_atoms = 4
box_dim = 10
timestep = 0.0004
temp = 100
fill=0

np.random.seed(1)
#init_vel = np.array([[0,0],[0,1],[1,0],[0.5,0.5]])
#init_pos = np.array([[5,5],[5 - 2**(1/6),5],[5,6],[4,4]])
x, num_atoms = fcc_lattice(num_atoms, box_dim, dim, fill)
init_vel, sigma= init_velocity(num_atoms, temp, dim)
init_pos = x
x, v, T, U, r = simulate(init_pos, init_vel, num_tsteps, timestep, box_dim, num_atoms, dim, temp)

'''
SHOWING VELOCITY DISTRIBUTION IS RANDOM:
mu, sigma = 0, 50
init_vel=np.random.normal(mu, sigma, size=(1000, dim))
plt.xlabel('init velocity')
s=init_vel.flatten()
count, bins, ignored = plt.hist(s, 30, density=True)
plt.plot(bins,1/(sigma*np.sqrt(2 * np.pi))*np.exp( - (bins - mu)**2 / (2 * sigma**2) ),linewidth=2, color='r')
plt.show()'''

# Animation and plotting stuff
if(dim == 2): # The animation only works for 2 dimensions.
    frames = num_tsteps
    fig, ax = plt.subplots()
    for k in range(len(x[0])):
        ax.plot(x[0,k,0], x[0,k,1], 'r.')
    ax.set_xlim(0, box_dim)
    ax.set_ylim(0, box_dim)

    def animate(i):
        ax.clear()
        ax.set_xlim(0, box_dim)
        ax.set_ylim(0, box_dim)
        for k in range(len(x[0])):
            ax.plot(x[i,k,0], x[i,k,1], 'r.')
    
    anim = matplotlib.animation.FuncAnimation(fig, animate, frames=frames, interval=1)
    anim
    #anim.save('Figures-animations/test2.gif', writer='imagemagick', fps=30)

plt.figure(2)
t = timestep * np.arange(0,num_tsteps)
plt.plot(t, U, label = 'potential energy')
plt.plot(t, T, label = 'kinetic energy')
plt.plot(t, T + U, label = 'total energy')
plt.xlabel(r'Time (in units of $\sqrt{\frac{m\sigma^2}{\epsilon}}$)')
plt.ylabel(r'Energy (in units of $\epsilon$)')
plt.legend()

#plt.figure(3)
#plt.plot(t, r)
plt.show()

print('specific heat computation:'')
print(skeleton_ours.specific_heat(T,num_atoms))