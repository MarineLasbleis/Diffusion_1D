

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from diffusion_1D import *


Nz = 50#number of space steps
diffusivity = 1.
Nt = 20 #number of time steps printed
timeMax = 20.

zmax = 0.
zmin = -10. #depth after wich we don't calculate

BC_0 = {'type':'Dir', 'T0':1., 'position':0}
BC_1 = {'type':'Dir', 'T0':0., 'position':-1}

# Initial conditions
z = np.linspace(zmin, zmax, Nz)
quantity =  0.*np.ones(Nz)
#quantity[0] = 0.
#quantity[-1] = 1.
quantity[BC_0['position']] = BC_0['T0']
quantity[BC_1['position']] = BC_1['T0']

plt.plot(z, quantity)

quantity_stock = np.zeros((Nz, Nt+1))
quantity_stock[:,0] = quantity


time = 0.
compteur_print = 1.

dx = np.amin(np.diff(z))
dt = 0.5*dx**2./diffusivity

Niterations = (timeMax - time)/dt
Ndiffplot = Niterations/Nt
print "test", Ndiffplot

nn = 0
while time < timeMax:

    nn += 1
    
    TCMB = np.sin(time*2.*np.pi/9.)+2* np.sin(time*2.*np.pi/4.)+0.5*np.sin(time*2.*np.pi/0.5)+3*np.sin(time*2.*np.pi/20.)
    BC_0['T0']=TCMB
    # compute new quantity (temperature/composition) 
    new_quantity = Solve_diffusion_cartesian(z, quantity, dx, dt, BC_0, BC_1, diffusivity)
    time = time+dt
    if int(nn%Ndiffplot) == 0 :
        print "plot number: ", int(nn/Ndiffplot)
        print "time: ", time, ". Iteration number: ", nn
        plt.plot(z, new_quantity)
        quantity_stock[:,int(nn/Ndiffplot)] = new_quantity
    quantity = new_quantity



fig, ax = plt.subplots(2)
time = np.linspace(0.,timeMax, Nt+1)
ax[0].plot(time, np.sin(time*2.*np.pi/9.)+2* np.sin(time*2.*np.pi/4.)+
             0.5*np.sin(time*2.*np.pi/0.5)+3*np.sin(time*2.*np.pi/20.))
ax[1].matshow(quantity_stock)



plt.show()


