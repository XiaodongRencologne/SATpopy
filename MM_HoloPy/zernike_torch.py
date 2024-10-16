
# coding: utf-8

# In[3]:


import numpy as np
import torch as T;
from math import factorial as fact
DEVICE=T.device('cuda' if T.cuda.is_available() else 'cpu');


# In[4]:


def ev(n, m, x, y):
    """Evaluate zernike polynomial (n,l)
    
    :param xy: The array over which to evaluate, of form [ [x,y] ]
    
    Example: evalZern(2,0, numpy.moveaxis(numpy.mgrid[-1:1:64j, -1:1:64j], 0, -1) )
    
    """
    M=abs(m)
    Nradterms=(n-M)//2+1;
    radpowers=n-2*np.arange(Nradterms); # rad powers
    radcoeffs= np.array(list(map(lambda s: (-1)**s * fact(n-s) / ( fact(s) * fact ( (n+M)/2 -s ) * fact( (n-M)/2 -s )),
                                    np.arange(Nradterms))))
    
    r=np.hypot(x,y);
    phi=np.arctan2(y,x);
    
    v=np.zeros_like(r);
    for i in range(Nradterms):
        v+=radcoeffs[i]*r**(radpowers[i]);
    if m>0:
        v*=np.cos(M*phi);
    elif m<0:
        v*=np.sin(M*phi);

    return v

def N(nmax):
    "Number of polynomials up to and including order nmax"
    return nmax*(nmax+3)//2 +1


# In[5]:


def mkCFn(N_order,x,y,dtype='numpy',device=T.device('cpu')):
    x=np.array(x);y=np.array(y);
    zz=np.array([]);
    for n in range(N_order+1):
        for m in range(-n,n+1,2):
            zz=np.append(zz,ev(n,m,x,y));
    N_vector=N(N_order);
    zz=zz.reshape(N_vector,-1);
    if dtype.lower()=='numpy':
        pass;
    elif dtype.lower()=='torch':
        zz=T.tensor(zz,dtype=T.float64).to(device);
    else:
        print('wrong dtype input!')
    
    
    def error(coeffs): 
        return (coeffs.reshape(N_vector,1)*zz).sum(0);
    error.parnames=["z%i"%i for i in range(zz.shape[0])]
    return error;
        
        

