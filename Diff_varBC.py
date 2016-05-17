

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from diffusion_1D import *


Nz = 50#number of space steps
diffusivity = 1.
Nt = 100 #number of time steps printed
timeMax = 10.

Nplots = 100 

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

quantity_stock = np.zeros((Nz, Nt))
quantity_stock[:,0] = quantity


time = 0.
compteur_print = 1.

dx = np.amin(np.diff(z))
dt = 0.5*dx**2./diffusivity

Niterations = (timeMax - time)/dt
Ndiffplot = Niterations/Nplots
print "test", Ndiffplot

nn = 0
while time < timeMax:

    nn += 1
    print "iteration number: ", nn
    print "time: ", time
    
    TCMB = np.sin(time*2.*np.pi/9.)
    BC_0['T0']=TCMB
    # compute new quantity (temperature/composition) 
    new_quantity = Solve_diffusion_cartesian(z, quantity, dx, dt, BC_0, BC_1, diffusivity)
    time = time+dt
    print nn%Ndiffplot
    if int(nn%Ndiffplot) == 0 :
        plt.plot(z, new_quantity)
    quantity = new_quantity






print dt, dx



plt.show()


