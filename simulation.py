import skeleton_ours
from skeleton_ours import simulate
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation
import time

num_tsteps = 20000
num_atoms = 4
box_dim = 10
timestep = 0.0004

np.random.seed(1)
init_vel = np.array([[0,0],[0,1],[1,0],[0.5,0.5]])
init_pos = np.array([[5,5],[5 - 2**(1/6),5],[5,6],[4,4]])
x, v, T, U, r = simulate(init_pos, init_vel, num_tsteps, timestep, box_dim, num_atoms)

"""
N = 1000
UT = np.zeros(N)
for i in range(N):
    x, v, T, U = simulate(init_pos, init_vel, num_tsteps, timestep, box_dim)
    UT[i] = U[0]

plt.plot(UT)
"""

# Animation and plotting stuff
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

# print(T)
# print(U)

plt.figure(2)
t = timestep * np.arange(0,num_tsteps)
plt.plot(t, U, label = 'potential energy')
plt.plot(t, T, label = 'kinetic energy')
plt.plot(t, T + U, label = 'total energy')
plt.xlabel(r'Time (in units of $\sqrt{\frac{m\sigma^2}{\epsilon}}$)')
plt.ylabel(r'Energy (in units of $\epsilon$)')
plt.legend()

plt.figure(3)
plt.plot(t, r)
plt.show()