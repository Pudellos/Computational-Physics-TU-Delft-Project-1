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

C=np.zeros((4, len_p))
D=np.zeros((len_p,dim))
D[0,0]=1
result=10

for i in range (1,len_p):
    C[0,i]=np.any([np.all(D-[D[i-1,0],D[i-1,1]+1]==0, axis=1)])
    C[1,i]=np.any([np.all(D-[D[i-1,0],D[i-1,1]-1]==0, axis=1)])
    C[2,i]=np.any([np.all(D-[D[i-1,0]+1,D[i-1,1]]==0, axis=1)])
    C[3,i]=np.any([np.all(D-[D[i-1,0]-1,D[i-1,1]]==0, axis=1)])
    C[:,i]=(1-C[:,i])
    A=np.random.rand(4)
    C[:,i]=C[:,i]*A
    result = np.argmax(C[:,i])
    if result==0:
        D[i,1]=D[i-1,1]+1
        D[i,0]=D[i-1,0]
    if result==1:
        D[i,1]=D[i-1,1]-1
        D[i,0]=D[i-1,0]
    if result==2:
        D[i,0]=D[i-1,0]+1
        D[i,1]=D[i-1,1]
    if result==3:
        D[i,0]=D[i-1,0]-1
        D[i,1]=D[i-1,1]
        

plt.figure(1)
plt.plot(D[:,0], D[:,1], label = 'Snake')