# Time-stamp: <2015-07-23 11:41:42 marine>


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from diffusion_1D import *


n=100 #number of space steps for each segments
#the time step will be obtained by specifiying CFL=1
kappa=0.01

Nt=10 #number of timesteps that will be printed



xAjouts=[1.,2.,3.]
timeAjouts=[0.,0.1,0.2]
tempAjouts=[1.,2.,3.]
Ajouts=dict(x=xAjouts, T=tempAjouts, time=timeAjouts)


#Initial conditions : 
x=np.linspace(0,Ajouts["x"][0],n)
T=np.ones(n)*Ajouts["T"][0]
T[0]=0.
T[-1]=Ajouts["T"][1]
plt.plot(x,T,'o')



def Solve_diffusion_cartesian(x,T,dx,dt,BC1,BC2,BCtype,kappa):

    a,b,c,d = CrankNicholson_dr2_order2(T,dx,kappa)
    b = b+1/dt
    d = d+T/dt
    d = BoundaryConditions(T,a,b,c,d,BC1,'bottom',typeBC=BCtype)
    d = BoundaryConditions(T,a,b,c,d,BC2,'top',typeBC=BCtype)
    T = TDMAsolver(a, b, c, d)
    
    return T

## Start


time=0.
compteur=1
screenshot=Ajouts["time"][1]/Nt


dx=np.amin(np.diff(x))
dt=dx**2./kappa

print dt, dx

while time<Ajouts["time"][1]:

    if time+dt>Ajouts["time"][1]*compteur/Nt:
        print time+dt, dt, Ajouts["time"][1]*compteur/Nt
        dt=Ajouts["time"][1]*compteur/Nt-time
        if dt<1e-12: dt=1e-12  
        time=Ajouts["time"][1]*compteur/Nt
        picture=1
        compteur=compteur+1
    else:
        picture=0
        dt=dx**2./kappa

    time=time+dt
    T=Solve_diffusion_cartesian(x,T,dx,dt,0.,0.,"Neumann",kappa,1.)

    if picture:
        print "picture", time, "CFL", kappa*dt/dx**2
        plt.plot(x,T)

x_2=np.linspace(Ajouts["x"][0]+dx/2.,Ajouts["x"][1],n)
x=np.concatenate((x,x_2),axis=1)
T=np.concatenate((T,Ajouts["T"][1]*np.ones(n)), axis=1)
print time
plt.show()
plt.plot(x,T,'o')

compteur=1
dx=np.amin(np.diff(x))

while time<Ajouts["time"][2]:

    if time+dt>Ajouts["time"][2]*compteur/Nt:
        print time+dt, dt, Ajouts["time"][2]*compteur/Nt
        dt=Ajouts["time"][2]*compteur/Nt-time    
        time=Ajouts["time"][2]*compteur/Nt
        picture=1
        compteur=compteur+1
    else:
        picture=0
        dt=dx**2./kappa

    time=time+dt
    T=Solve_diffusion_cartesian(x,T,dx,dt,0.,Ajouts["T"][2],"Dir",kappa)

    if picture:
        print "picture", time
        plt.plot(x,T),"CFL", kappa*dt/dx**2
   
    

plt.show()
