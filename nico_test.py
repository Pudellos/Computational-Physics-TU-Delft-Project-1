import numpy as np
import matplotlib.pyplot as plt
from numpy import random

def force(i,x,size):#returns force vector (d dimensions) acting on particle i given the positions of all particles (x)
    e = 1#119.8*1.38*(10**-23) #epsilon
    S = 1#3.405*(10**-10) #sigma
    xr = np.concatenate((x[0:i,:],x[i+1:len(x),:])) #remove particle i from matrix x to aviod self-interaction
    deltax=x[i,:]-xr
    deltax=np.where(deltax > (-size/2),deltax, deltax-size) #distances smaller than -size/2 get - size term
    deltax=np.where(deltax < (size/2),deltax, deltax+size)  #distances larger than size/2 get + size term
    r = np.power(np.sum(np.power(deltax,2),1),1/2) #calculate distance between i and all particles
    F = -np.dot(np.transpose(deltax),(48*e*(np.power(S,12)/np.power(r,14)-0.5*np.power(S,6)/np.power(r,8)))) #sum of all Lenard-Jones vectors
    return F


def simulate(init_pos, init_vel, num_tsteps, timestep, box_dim, size):
    x0 = init_pos #assign random initial particle positions
    N = len(x0) #particle number
    x = np.zeros((N,box_dim))+x0
    v0 = init_vel #assign random initial particle velocities
    v = np.zeros((N,box_dim))+v0
    F = np.zeros((N,box_dim))
    m = 1#6.6335209*(10**-26) #mass
    for i in range(num_tsteps): #time evolution
        x=x+v*timestep #update position
        v=v+F*timestep/m #update velocity
        for k in range(N): #calculate force acting on every particle
            F[k]=force(k,x,size)
            for a in range(box_dim): #apply periodic boundary conditions (in d dimensions) for all particles
                if x[k,a] > size:
                    x[k,a]=x[k,a]-size
                if x[k,a] < 0:
                    x[k,a]=x[k,a]+size
    return x,v

def energy(x,v): #returns total energy in the system (Lenard-Jones potential + kinetic)
    e = 1#119.8*1.38*(10**-23) #epsilon
    S = 1#3.405*(10**-10) #sigma
    m = 1#6.6335209*(10**-26) #mass
    ULJ=0
    EV=0
    for i in range(len(x)):
        xr = np.concatenate((x[0:i,:],x[i+1:len(x),:])) #remove particle i from matrix x to aviod self-interaction
        r = np.power(np.sum(np.power(xr-x[i,:],2),1),1/2) #calculate distance between i and all particles
        ULJ = ULJ + sum(4*e*(np.power(S,12)/np.power(r,12)-np.power(S,6)/np.power(r,6))) #sum of all L-J potentials for partical i
    EV=0.5*m*sum(np.sum(np.power(v,2),1)) #sum of the the kinetic energy of all particles
    return ULJ, EV