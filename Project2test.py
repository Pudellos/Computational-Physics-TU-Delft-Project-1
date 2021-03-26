import itertools
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation
import time
from scipy.optimize import curve_fit
import random
dim=2
num_p=20
len_p=20

Allow_D=np.zeros((4, len_p))
Snake=np.zeros((len_p,dim))
Snake[0,0]=1
Chosen_D=10

for i in range (1,len_p):
    Allow_D[0,i]=np.any([np.all(Snake-[Snake[i-1,0],Snake[i-1,1]+1]==0, axis=1)])    #determines which direction one can choose from
    Allow_D[1,i]=np.any([np.all(Snake-[Snake[i-1,0],Snake[i-1,1]-1]==0, axis=1)])    
    Allow_D[2,i]=np.any([np.all(Snake-[Snake[i-1,0]+1,Snake[i-1,1]]==0, axis=1)])
    Allow_D[3,i]=np.any([np.all(Snake-[Snake[i-1,0]-1,Snake[i-1,1]]==0, axis=1)])
    Allow_D[:,i]=(1-Allow_D[:,i])
    A=np.random.rand(4) #random direction
    Allow_D[:,i]=Allow_D[:,i]*A
    Chosen_D = np.argmax(Allow_D[:,i])   #chooses from allowed direction
    if Chosen_D==0:
        Snake[i,1]=Snake[i-1,1]+1
        Snake[i,0]=Snake[i-1,0]
    if Chosen_D==1:
        Snake[i,1]=Snake[i-1,1]-1
        Snake[i,0]=Snake[i-1,0]
    if Chosen_D==2:
        Snake[i,0]=Snake[i-1,0]+1
        Snake[i,1]=Snake[i-1,1]
    if Chosen_D==3:
        Snake[i,0]=Snake[i-1,0]-1
        Snake[i,1]=Snake[i-1,1]
        

plt.figure(1)
plt.plot(Snake[:,0], Snake[:,1], label = 'Snake')