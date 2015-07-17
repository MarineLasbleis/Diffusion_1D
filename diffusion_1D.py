### Diffusion_1D
# Time-stamp: <2015-07-17 10:35:47 marine>

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
    '''
    Solve the equation ax_{i-1}+bx_{i}+cx_{i+1}=d_{i} for x
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    https://gist.github.com/ofan666/1875903

    requisite:
    a_1=0
    c_end=0
    '''
    nf = len(a)     # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d))     # copy the array
    for it in xrange(1, nf):
        mc = ac[it]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
 
    xc = ac
    xc[-1] = dc[-1]/bc[-1]
 
    for il in xrange(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]
 
    del bc, cc, dc  # delete variables from memory
 
    return xc



def CrankNicholson_dr2_order2(T,dr,k):
    '''
    Compute the terms of the CrankNicholson scheme for the term d^2T/dr^2

    
    T : temperature profile
    dr and k are constants:
    a=k/2dr^2,    b=-k/dr^2,    c=k/2dr^2,     d=-k/2dr^2*(T(i-1)-2T(i)+T(i+1))
    '''

    nf = len(T)

    a=-k/2./dr**2*np.ones(nf)
    b=k/dr**2*np.ones(nf)
    c=-k/2./dr**2*np.ones(nf)

    d=np.zeros(nf)
    for i in xrange(1,nf-1):
        d[i]=k/2./dr**2*(T[i-1]-2.*T[i]+T[i+1])
    d[0], d[-1]=0.,0.
        
    return a, b, c, d


def BoundaryConditions(T,a,b,c,d,BC,pos,type='Dir'):
    '''
    Modify the a,b,c,d given the type of boundary conditions (Dirichlet by default)
    at the position specified by pos (top or bottom)

    BC is the value of the boundary conditions (a temperature in case of Dirichlet conditions, a flux in case of Neumann conditions)
    '''

    if type=='Dir':
        ## if pos=="top":
        ##     d[-2] = d[-2] -c[-2]*BC
        ##     a[-1], b[-1], c[-1], d[-1] = 0.,1.,0., BC
        ## if pos=="bottom":
        ##     d[1] = d[1] -c[1]*BC
        ##     a[0], b[0], c[0], d[0] = 0., 1., 0., BC

        if pos=="top":
            d[-1] = d[-1] -c[-1]*BC        
        if pos=="bottom":
            d[0] = d[0] -c[0]*BC
          
                
    return a,b,c,d



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




## Solver : test 1 : temperature in a dyke of thickness Lx=1

if Test1:
    ### Temperature profile

    n=1000 #size of the system
    x=np.linspace(0.5, 2.5, n)
    dr=x[1]-x[0]
    T0=0.5 #initial conditions
    k=1.e-6 #diffusion
    tc=1/k/np.pi**2
    nt=1000
    time=np.linspace(0,0.1*tc,nt)
    dt=time[1]-time[0]
    nplots=200

    ## Initial conditions
    T=T0*sum_truncated_init(x,10000,Lx=1.)
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
        plt.plot(x,T,'o')
        plt.plot(x,T0*sum_truncated_sol(x,10000,time[it],k),'.')
        print "time :", it, time[it]

    print 'CFL', k*dt/dr**2
    plt.show()
