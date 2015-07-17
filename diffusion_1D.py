### Diffusion_1D
# Time-stamp: <2015-07-17 10:47:17 marine>

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






