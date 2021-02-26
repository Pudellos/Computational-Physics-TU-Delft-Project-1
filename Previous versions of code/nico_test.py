import numpy as np
from numpy import random

def force(i,x,size):#returns force vector (d dimensions) acting on particle i given the positions of all particles (x)
    xr = np.concatenate((x[0:i,:],x[i+1:len(x),:])) #remove particle i from matrix x to aviod self-interaction
    deltax=x[i,:]-xr
    deltax=np.where(deltax < (-size/2),deltax+size, deltax) #distances smaller than -size/2 get + size term
    deltax=np.where(deltax > (size/2),deltax-size, deltax)  #distances larger than size/2 get - size term
    r = np.power(np.sum(np.power(deltax,2),1),1/2) #calculate distance between i and all particles
    F = np.dot(np.transpose(deltax),(48*(1/np.power(r,14)-0.5*(1/np.power(r,8))))) #sum of all nondimensional Lenard-Jones vectors       
    return F, deltax

def simulate(init_pos, init_vel, num_tsteps, timestep, box_dim, size):
    '''This is non-dimensional or in other words the input, init_pos is in units of sigma
    the Lenard Jones potential is in units of epsilon and time is in units of sigma*sqrt(m/e)
    sigma = 3.405*(10**-10) meter
    mass = 6.6335209*(10**-26) kilogram    
    epsilon = 119.8*1.38*(10**-23) K/kb'''

    x0 = init_pos #assign initial particle positions
    N = len(x0) #particle number
    x = np.zeros((N,box_dim))+x0
    v0 = init_vel #assign random initial particle velocities
    v = np.zeros((N,box_dim))+v0
    F = np.zeros((N,box_dim))
    dx= np.zeros((num_tsteps,1))
    ULJ= np.zeros((num_tsteps,1))
    EV= np.zeros((num_tsteps,1))
    for i in range(num_tsteps): #time evolution
        dx[i]= np.power(np.sum(np.power(x[0,:]-x[1,:],2)),1/2)
        U, V = energy(x,v)
        ULJ[i]=U
        EV[i]= V
        x=x+v*timestep #update position
        v=v+F*timestep #update velocity
        x=np.where(x>size,x-size,x) #apply periodic BC
        x=np.where(x<0,x+size,x) #apply periodic BC
        for k in range(N): #calculate force acting on every particle
            F[k], deltax = force(k,x,size)
    return x,v,dx,deltax,ULJ,EV

def energy(x,v): #returns total energy in the system (Lenard-Jones potential + kinetic)
    ULJ=0
    EV=0
    for i in range(len(x)):
        xr = np.concatenate((x[0:i,:],x[i+1:len(x),:])) #remove particle i from matrix x to aviod self-interaction
        r = np.power(np.sum(np.power(xr-x[i,:],2),1),1/2) #calculate distance between i and all particles
        ULJ = ULJ + sum(4*(1/np.power(r,12)-1/np.power(r,6))) #sum of all L-J potentials for partical i
    EV=0.5*sum(np.sum(np.power(v,2),1)) #sum of the the kinetic energy of all particles
    return ULJ, EV