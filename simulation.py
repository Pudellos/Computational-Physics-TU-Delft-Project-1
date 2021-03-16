import skeleton_ours
from skeleton_ours import simulate, fcc_lattice, init_velocity, normal_autocorr, autocorrelation, func, specific_heat, mean_squared_displacement
import itertools
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation
import time
from scipy.optimize import curve_fit

dim = 3
num_tsteps = 4000
num_atoms = 16
box_dim = 10
timestep = 0.004
temp = 100
fill = 0

x, num_atoms = fcc_lattice(num_atoms, box_dim, dim, fill)
init_vel, sigma = init_velocity(num_atoms, temp, dim)
init_pos = x
x, v, K, U = simulate(init_pos, init_vel, num_tsteps, timestep, box_dim, num_atoms, dim, temp)

if(False):
    MSD, AMSD, D, TSD1 = mean_squared_displacement(x, box_dim)
    print("done")
    AC = autocorrelation(AMSD)
    print(" AC done")

    plt.subplot(121)
    plt.plot(AMSD)
    plt.ylabel(r'$\langle\Delta^2x(t)\rangle$')
    plt.xlabel(r'Time $t$')

    plt.subplot(122)
    plt.plot(AC)
    plt.ylabel(r'$\chi(\tau)$')
    plt.xlabel(r'Time lag $\tau$')
    plt.show()

# Code to calculate heat capacity
if(True):
    #np.random.seed(1)
    Cv = np.zeros(100)
    for i in range(100):
        x, num_atoms = fcc_lattice(num_atoms, box_dim, dim, fill)
        init_vel, sigma = init_velocity(num_atoms, temp, dim)
        init_pos = x
        x, v, K, U = simulate(init_pos, init_vel, num_tsteps, timestep, box_dim, num_atoms, dim, temp)

        plt.plot(K)
        Cv[i] = specific_heat(K, num_atoms)
        print('Specific heat Cv:', Cv[i])

    print(np.mean(Cv))
    print(np.std(Cv))
    plt.show()
    
"""
plt.figure(2)
t = timestep * np.arange(0,len(U))
plt.plot(t, U, label = 'potential energy')
plt.plot(t, K, label = 'kinetic energy')
plt.plot(t, K + U, label = 'total energy')
plt.xlabel(r'Time (in units of $\sqrt{\frac{m\sigma^2}{\epsilon}}$)')
plt.ylabel(r'Energy (in units of $\epsilon$)')
plt.legend()
plt.show()
"""

## Testing of autocorrelation function ()
if(False):
    correlation_time = 50
    t_max = 20000
    testcor = normal_autocorr(0, 1, correlation_time, t_max)
    
    Xi = autocorrelation(testcor)
    t = np.arange(0,t_max)
    lastfit = 100
    t_sub = np.arange(0,lastfit)
    popt, pcov = curve_fit(func, t_sub, Xi[0:lastfit], bounds=(0, [10.]))

    #Xi = np.mean(XIs, axis=0)
    #plt.plot(t_sub, np.exp(-popt[0] * t_sub))
    #plt.plot(t_sub, np.exp(-t_sub / correlation_time))
    #plt.plot(t, Xi)
    print("The fitted correlation time is ", 1/popt[0],"the actual correlation time is ", correlation_time)

    plt.figure(1)
    plt.plot(testcor)

    plt.figure(2)
    plt.plot(t, Xi)
    
    plt.show()

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