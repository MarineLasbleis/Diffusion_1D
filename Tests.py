import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import diffusion_1D



Test1=1





def sum_truncated_init(x,pmax,Lx=1.):
    '''
    Compute the sum for the initial condition (see Dubuffet course p22), troncated at index pmax
    '''
    somme=np.zeros(len(x))
    for p in xrange(0,pmax):
        somme=somme+np.sin((2.*p+1.)*np.pi*(x+1.)/Lx)/(2.*p+1.)
    return 1.+ 4./np.pi*somme

def sum_truncated_sol(x,pmax,t,kappa,Lx=1.):
    '''
    Compute the sum for the initial condition (see Dubuffet course p22), troncated at index pmax
    '''
    somme=np.zeros(len(x))
    for p in xrange(0,pmax):
        somme=somme+np.sin((2.*p+1.)*np.pi*(x+1.)/Lx)/(2.*p+1.)*np.exp(-1*(2.*p+1.)**2*np.pi**2*kappa*t/Lx**2)
    return 1.+ 4./np.pi*somme




if Test1:
    
    n=1000 #size of the system
    x=np.linspace(0.5, 2.5, n)
    dr=x[1]-x[0]
    T0=0.5 #initial conditions
    k=1.e-6 #diffusion
    tc=1/k/np.pi**2
    nt=3000
    time=np.linspace(0,0.1*tc,nt)
    dt=time[1]-time[0]
    nplots=600
    Lx=1.

    print "Test : diffusion of a rectangular function of thickness", Lx
    print "The code solve the diffusion equation with Dirichlet BC, and plot also the analytical solution"
    print nt/nplots-1, ' times steps are plotted on the graph, with ', dt*nplots, 'seconds btwn each plots'

    ## Initial conditions
    T=T0*sum_truncated_init(x,10000,Lx)
    plt.plot(x,T)

    for it in xrange(1,nt):

        a,b,c,d = CrankNicholson_dr2_order2(T,dr,k)
        b = b+1/dt
        d = d+T/dt

        a,b,c,d = BoundaryConditions(T,a,b,c,d,0.,"top",type='Dir')
        a,b,c,d = BoundaryConditions(T,a,b,c,d,0.,"bottom",type='Dir')

        #a[0],b[0],c[0],d[0] = 0., -1.,1.,0.
        #a[-1],b[-1],c[-1],d[-1]=-1.,1.,0.,0.

        xc = TDMAsolver(a, b, c, d)   
        T=xc

        if (it % nplots)== 0 :
            plt.plot(x[0:-1:20],T[0:-1:20],'o')
            plt.plot(x,T0*sum_truncated_sol(x,10000,time[it],k))
            print "time :", it, time[it]

    print 'CFL', k*dt/dr**2
    plt.show()
