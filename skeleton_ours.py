import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation
import time
import itertools

plt.rcParams['figure.dpi']=100
plt.rcParams["animation.html"] = "jshtml"
plt.rcParams['animation.convert_path'] = "F:/Python"

### Constants ###
k_B = 1.380649e-23  #J K^-1
eps = 119.8 * k_B   # J
sigma = 3.405       # Angstrom
m = 6.6e-26         # kg

def simulate(init_pos, init_vel, num_tsteps, timestep, box_dim, num_atoms, dim, temp):
    """
    Molecular dynamics simulation using the Euler or Verlet's algorithms
    to integrate the equations of motion. Calculates energies and other
    observables at each timestep.

    Parameters
    ----------
    init_pos : np.ndarray
        The initial positions of the atoms in Cartesian space
    init_vel : np.ndarray
        The initial velocities of the atoms in Cartesian space
    num_tsteps : int
        The total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    Any quantities or observables that you wish to study.
    x : all positions
    v : all velocities
    E : energy at every time step
    p : all momenta
    """
    h = timestep # Shorther way to deal with timestep :)
    x = np.ndarray(shape = (num_tsteps, num_atoms, dim)) # Stores velocities for every timestep
    v = np.ndarray(shape = (num_tsteps, num_atoms, dim)) # Stores positions for every timestep
    F = np.ndarray(shape = (num_tsteps, num_atoms, dim)) # Stores force for every timestep

    T = np.zeros(num_tsteps) # Kinetic energy
    U = np.zeros(num_tsteps) # Potential energy    
    
    # Initialize position, velocity and force
    x[0] = init_pos
    v[0] = init_vel

    # Calculate the energies and force for thwe first timestep
    rel_pos, rel_dist = atomic_distances(x[0], box_dim)
    U[0] = potential_energy(rel_dist)
    T[0] = kinetic_energy(v[0])
    F[0] = lj_force(rel_pos, rel_dist, dim)

    rescale = 0
    eq_time = 0
    for i in range(1, num_tsteps):
        
        # Verlet method to calculate new position and velocity
        x[i] = x[i - 1] + v[i - 1] * h + (h**2 / 2) * F[i - 1]
        rel_pos, rel_dist = atomic_distances(x[i], box_dim)
        F[i] = lj_force(rel_pos, rel_dist, dim)
        v[i] = v[i - 1] + (h / 2) * ( F[i] + F[i - 1] )
        
        #Euler method to calculate new position and velocity
        #x[i] = x[i - 1] + v[i - 1] * h
        #v[i] = v[i - 1] + F * h

        # Periodic boundary conditions
        for d in range(dim):
            x[i, :, d] = np.where( x[i, :, d] < box_dim  , x[i, :, d], x[i, :, d] % box_dim)
            x[i, :, d] = np.where( x[i, :, d] > 0        , x[i, :, d], x[i, :, d] % box_dim)

        # Calculate the energies
        T[i] = kinetic_energy(v[i])
        U[i] = potential_energy(rel_dist)
    
        # After c steps check if the target temperature is reached, otherwise rescale velocities
        c = 20         
        if((i % c) == 0 and rescale == 0): 
                E_kin = (num_atoms - 1) * (3 / 2) * (k_B / eps) * temp
                E_avg = np.mean(T[i-c:i])
                L = np.sqrt(E_kin/ E_avg) 
                v[i] = L * v[i]
                
                if(np.abs(E_kin - E_avg) < 0.05):
                    #print("Temperature is: ", E_avg / ((num_atoms - 1) * (3 / 2) * (k_B / eps)) )
                    rescale = 1                    
                    eq_time = i

    return (x[eq_time:], v[eq_time:], T[eq_time:], U[eq_time:])

def atomic_distances(pos, box_dim):
    """
    Calculates relative positions and distances between particles.

    parameters
    ----------
    pos : np.ndarray
        The positions of the particles in cartesian space
    box_dim : float
        The dimension of the simulation box

    returns
    -------
    rel_pos : np.ndarray
        Relative positions of particles
    rel_dist : np.ndarray
        The distance between particles
    """
    rel_pos = pos[:, np.newaxis] - pos

    rel_pos = np.where(rel_pos < (-box_dim / 2), rel_pos + box_dim, rel_pos) #distances smaller than -size/2 get + size term
    rel_pos = np.where(rel_pos > ( box_dim / 2), rel_pos - box_dim, rel_pos)  #distances larger than size/2 get - size ter

    rel_dist = np.linalg.norm(rel_pos, axis=2)    
    return (rel_pos, rel_dist)


def lj_force(rel_pos, rel_dist, dim):
    """
    Calculates the net forces on each atom.

    Parameters
    ----------
    rel_pos : np.ndarray
        Relative particle positions as obtained from atomic_distances
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    np.ndarray
        The net force acting on particle i due to all other particles
    """
    r = np.array([rel_dist]*dim).transpose()
    F = rel_pos * ( 48 * np.power(r, -14, where= r!=0) - 24 * np.power(r, -8, where= r!=0) )
    F = np.sum(F, axis=1)
    return F


def fcc_lattice(num_atoms, box_dim, dim, fill):
    """
    Initializes a system of atoms on an fcc lattice.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system
    box_dim : float
        The dimension of the simulation box
    dim : int
        The number of dimensions we are looking at
    fill: int
        If fill=0 the fcc lattice will only use num_atoms particles,
        If fill is not 0 the fcc lattice will be completely filled
    Returns
    -------
    x : np.ndarray
        Array of particle coordinates
    len(x) : int
        Number of particles in the fcc lattice
    """
    ppuc = 2 * dim - 2 #particles per unit cell (4 in 3d and 2 in 2d)
    UC = num_atoms / ppuc #total number of unit cells needed
    AUC = np.ceil(UC**(1 / dim)) #number of unit cells per axis
    lc = box_dim / AUC #lattice constant
    L = np.arange(AUC)
    comb = [p for p in itertools.product(L, repeat=dim)] #all possible permutation to fill the whole lattice (all unit cells once)
    a = np.asarray(comb)
    if dim == 3:
        x = np.append(a,a+[0.5,0.5,0], axis=0) #add basisvectors to unit fcc cell (3d)
        x = np.append(x,a+[0.5,0,0.5], axis=0)
        x = np.append(x,a+[0,0.5,0.5], axis=0)
    if dim == 2:
        x = np.append(a,a+[0.5,0.5], axis=0) #add basisvectors to unit fcc cell (2d)
    x = lc * x
    if fill == 0:
        x=x[0:num_atoms]
    print("To fill the FCC lattice", len(x)-num_atoms, "particles were added for a total of", len(x),"particles.")
    return x, len(x)

def kinetic_energy(vel):
    """
    Computes the kinetic energy of an atomic system.

    Parameters
    ----------
    vel: np.ndarray
        Velocity of particle

    Returns
    -------
    float
        The total kinetic energy of the system.
    """
    T = np.linalg.norm(vel, axis=1)
    T = 1/2 * np.power(T, 2)
    return np.sum(T)

def potential_energy(rel_dist):
    """
    Computes the potential energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    float
        The total potential energy of the system.
    """
    r = np.copy(rel_dist)
    r[np.allclose(r, 0)] = np.inf
    U = 4 * (np.power(r, -12, where=~np.isclose(r,0)) - np.power(r, -6, where=~np.isclose(r,0)) )
    U = np.sum(U)    
    return U / 2

def init_velocity(num_atoms, temp, dim):
    """
    Initializes the system with Gaussian distributed velocities.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system.
    temp : float
        The (unitless) temperature of the system.

    Returns
    -------
    vel_vec : np.ndarray
        Array of particle velocities
    """    
    sigma = ((k_B/eps) * temp)**0.5
    return np.random.normal(0, sigma, (num_atoms, dim)), sigma

def specific_heat(T, num_atoms):
    '''
    Calculates specific heat under constant volume of the simulated gas.
    
    Parameters
    ----------
    T : np.ndarray
        kinetic energy
    num_atoms : int
        The total number of simulated atoms

    Returns
    -------
    Cv : float
        specific heat of the gas under constant volume condition'''
    
    mean = np.mean(T)
    fluct = np.var(T)
    Cv = 3 * num_atoms / (2 - 3 * num_atoms * (fluct / mean**2))
    return Cv

def mean_squared_displacement(x, box_dim):
    
    '''
    Calculates the mean squared displacement on each atom and the average. The first timestep is the reference position for the calculation.

    Parameters
    ----------
    x : np.ndarray
        Positions of all particles for all timesteps
    box_dim : np.ndarray
        Linear size of the domain
        
    Returns
    -------
    MSD : np.ndarray
        MSD(t) for all individual particles
    ASMD: np.ndarray
        Averaged MSD(t) over all particles
    '''
    TSD1 = x[1:len(x)] - x[0:len(x)-1]
    TSD1 = np.where(TSD1 >  0.9 * box_dim , TSD1 - box_dim, TSD1)
    TSD1 = np.where(TSD1 < -0.9 * box_dim , TSD1 + box_dim, TSD1)
    D = np.cumsum(TSD1, axis=0)
    MSD = np.sum(np.power(D,2), axis=2)
    AMSD = np.sum(MSD,1) / x.shape[1]
    return MSD, AMSD, D, TSD1

def func(x, b):
    return np.exp(- b * x)

def autocorrelation(A):    
    '''
    Calculates the autocorrelation function for a given dataset

    Parameters
    ----------
    A : np.ndarray
        The data to be correlated
        
    Returns
    -------
    Xi : array of the autocorrelation function(t)
    '''
    N = len(A)
    Xi = np.zeros(N - 1)

    for t in range(N - 1):
        #print("%s/%s" %(t,N))
        t1 = (N - t) * sum( A[0:(N - t)] * A[t:N]) #term 1
        t2 = sum(A[0:(N - t)]) * sum(A[t:N]) #term 2
        t3 = np.sqrt( (N - t) * sum(A[0:(N - t)]**2) - sum(A[0:(N - t)])**2 ) #term 3
        t4 = np.sqrt( (N - t) * sum(A[t:N]**2) - sum(A[t:N])**2 ) #term 4

        Xi[t] = ( t1 - t2 ) / ( t3 * t4 ) 
    return Xi

def normal_autocorr(mu, sigma, tau, N):
    """Generates an autocorrelated sequence of Gaussian random numbers.
    
    Each of the random numbers in the sequence of length `N` is distributed
    according to a Gaussian with mean `mu` and standard deviation `sigma` (just
    as in `numpy.random.normal`, with `loc=mu` and `scale=sigma`). Subsequent
    random numbers are correlated such that the autocorrelation function
    is on average `exp(-n/tau)` where `n` is the distance between random
    numbers in the sequence.
    
    This function implements the algorithm described in
    https://www.cmu.edu/biolphys/deserno/pdf/corr_gaussian_random.pdf
    
    Parameters
    ----------
    
    mu: float
        mean of each Gaussian random number
    sigma: float
        standard deviation of each Gaussian random number
    tau: float
        autocorrelation time
    N: int
        number of desired random numbers
    
    Returns:
    --------
    sequence: numpy array
        array of autocorrelated random numbers
    """
    f = np.exp(-1./tau)
    
    sequence = np.zeros(shape=(N,))
    
    sequence[0] = np.random.normal(0, 1)
    for i in range(1, N):
        sequence[i] = f * sequence[i-1] + np.sqrt(1 - f**2) * np.random.normal(0, 1)
    
    return mu + sigma * sequence