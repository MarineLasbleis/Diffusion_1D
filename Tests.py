import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys

from diffusion_1D import *


Test=sys.argv[1]





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




if Test=="1":
    
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

    print "==========================================================================="
    print "=== Test : diffusion of a rectangular function of thickness", Lx
    print "=== The code solve the diffusion equation with Dirichlet BC, and plot also the analytical solution"
    print "=== ", nt/nplots-1, ' times steps are plotted on the graph, with ', dt*nplots, 'seconds btwn each plots'

    ## Initial conditions
    T=T0*sum_truncated_init(x,10000,Lx)
    plt.plot(x,T)

    for it in xrange(1,nt):

        a,b,c,d = CrankNicholson_dr2_order2(T,dr,k)
        b = b+1/dt
        d = d+T/dt
        
        d = BoundaryConditions(T,a,b,c,d,0,'top',typeBC='Dir', dx=0.,BC2=0.)
        d = BoundaryConditions(T,a,b,c,d,0,'bottom',typeBC='Dir', dx=0.,BC2=0.)

        #a[0],b[0],c[0],d[0] = 0., -1.,1.,0.
        #a[-1],b[-1],c[-1],d[-1]=-1.,1.,0.,0.

        
        xc = TDMAsolver(a, b, c, d)   
        T=xc

        if (it % nplots)== 0 :
            plt.plot(x[0:-1:20],T[0:-1:20],'o')
            plt.plot(x,T0*sum_truncated_sol(x,10000,time[it],k))
            print "== time :", it, time[it]

    print '== CFL', k*dt/dr**2
    print '=== Figure  Test1'
    plt.show()
    print "==========================================================================="



elif Test=="2":

    ## solve the diffusion of T(x,0)=sin(pi *x/L)

    L=1.
    alpha=1.
    n=200
    nt=20000
    nplots=4000
    time=np.linspace(0,0.1,nt)
    dt=(time[1]-time[0])

    x1=np.linspace(0,L/2.,n)
    x2=np.linspace(L/2.,L,n)
    dx12=x1[1]-x1[0]
    x3=np.linspace(0,L,n)
    dx3=x3[1]-x3[0]
    
    
    T1=np.sin(x1*np.pi/L)
    T2=np.sin(x2*np.pi/L)
    T3=np.sin(x3*np.pi/L)

    f, ax=plt.subplots(1,3)
    ax[0].plot(x1,T1)
    ax[1].plot(x2,T2)
    ax[2].plot(x3,T3)

    ax[2].plot(x3,np.exp(-alpha*np.pi**2/L**2*time[0])*np.sin(np.pi*x3/L),'ob')
    ax[1].plot(x2,np.exp(-alpha*np.pi**2/L**2*time[0])*np.sin(np.pi*x2/L),'ob')
    ax[0].plot(x1,np.exp(-alpha*np.pi**2/L**2*time[0])*np.sin(np.pi*x1/L),'ob')
    
    
    for it in xrange(1,nt):

        a,b,c,d = CrankNicholson_dr2_order2(T1,dx12,alpha)
        b = b+1/dt
        d = d+T1/dt
        d = BoundaryConditions(T1,a,b,c,d,0.,'bottom',typeBC='Dir')
        d = BoundaryConditions(T1,a,b,c,d,0.,'top',typeBC='Neumann')

        xc = TDMAsolver(a, b, c, d)   
        T1=xc

        a,b,c,d = CrankNicholson_dr2_order2(T2,dx12,alpha)
        b = b+1/dt
        d = d+T2/dt
        d = BoundaryConditions(T2,a,b,c,d,0.,'bottom',typeBC='Neumann')
        d = BoundaryConditions(T2,a,b,c,d,0.,'top',typeBC='Dir')

        xc = TDMAsolver(a, b, c, d)   
        T2=xc
        
        a,b,c,d = CrankNicholson_dr2_order2(T3,dx3,alpha)
        b = b+1/dt
        d = d+T3/dt
        d = BoundaryConditions(T3,a,b,c,d,0.,'bottom',typeBC='Dir')
        d = BoundaryConditions(T3,a,b,c,d,0.,'top',typeBC='Dir')

        xc = TDMAsolver(a, b, c, d)   
        T3=xc

    

        if (it % nplots)== 0 :
            print it, T3[0]
            ax[0].plot(x1,T1)
            ax[1].plot(x2,T2)
            ax[2].plot(x3,T3)
            ax[2].plot(x3[::10],np.exp(-alpha*np.pi**2/L**2*time[it])*np.sin(np.pi*x3[::10]/L),'ob')
            ax[1].plot(x2[::10],np.exp(-alpha*np.pi**2/L**2*time[it])*np.sin(np.pi*x2[::10]/L),'ob')
            ax[0].plot(x1[::10],np.exp(-alpha*np.pi**2/L**2*time[it])*np.sin(np.pi*x1[::10]/L),'ob')

    print 'CFL - 1-2',  alpha*dt/dx12**2
    print 'CFL - 3',  alpha*dt/dx3**2
    
    plt.show()


else :
    print "You have to choose a valid test number."
