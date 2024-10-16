#%%
import numpy as np
import torch as T
import copy
DEVICE0=T.device('cpu')
# %%
'''
1. define a vector
'''
class vector():
    def __init__(self):
        self.x=np.array([])
        self.y=np.array([])
        self.z=np.array([])
    def to_Tensor(self,DEVICE='cpu'):
        '''
        convert numpy array into tensor and 
        load the data into cpu RAM or GPU RAM.
        '''
        Device=T.device(DEVICE)
        if type(self.x).__module__ == np.__name__:
            self.x=T.tensor(self.x,dtype=T.float64).to(Device)
            self.y=T.tensor(self.y,dtype=T.float64).to(Device)
            self.z=T.tensor(self.z,dtype=T.float64).to(Device)
        elif type(self.x).__module__ == T.__name__:
            self.x=self.x.to(Device)
            self.y=self.y.to(Device)
            self.z=self.z.to(Device)
    def to_numpy(self):
        if type(self.x).__module__ == np.__name__:
            pass
        elif type(self.x).__module__ == T.__name__:
            self.x=self.x.cpu().numpy()
            self.y=self.y.cpu().numpy()
            self.z=self.z.cpu().numpy()
        

        
'''
2.
'''
def To_co(theta,phi):
    r0=vector()
    theta0=vector()
    PHi0=vector()

    r0.x=np.sin(theta)*np.cos(phi)
    r0.y=np.sin(theta)*np.sin(phi)
    r0.z=np.cos(theta)
    
    theta0.x=np.cos(theta)*np.cos(phi)
    theta0.y=np.cos(theta)*np.sin(phi)
    theta0.z=-np.sin(theta)
    
    PHi0.x=-np.sin(phi)
    PHi0.y=np.cos(phi)
    PHi0.z=np.zeros(phi.size)

    co=sumvector(scalarproduct(np.cos(phi),theta0),scalarproduct(-np.sin(phi),PHi0))
    cx=sumvector(scalarproduct(np.sin(phi),theta0),scalarproduct(np.cos(phi),PHi0))
    crho=r0
    return co,cx,crho



'''
2. electromagnetic field vector
'''
class Fvector():
    def __init__(self,):
        self.x=np.zeros(0,dtype=np.complex128)
        self.y=np.zeros(0,dtype=np.complex128)
        self.z=np.zeros(0,dtype=np.complex128)
    def to_Tensor(self,DEVICE='cpu'):
        Device=T.device(DEVICE)
        if type(self.x).__module__ == np.__name__:
            self.x=T.tensor(self.x,dtype=T.cdouble).to(Device)
            self.y=T.tensor(self.y,dtype=T.cdouble).to(Device)
            self.z=T.tensor(self.z,dtype=T.cdouble).to(Device)
        elif type(self.x).__module__ == T.__name__:
            self.x=self.x.to(Device)
            self.y=self.y.to(Device)
            self.z=self.z.to(Device)
    def to_numpy(self):
        if type(self.x).__module__ == np.__name__:
            pass
        elif type(self.x).__module__ == T.__name__:
            self.x=self.x.cpu().numpy()
            self.y=self.y.cpu().numpy()
            self.z=self.z.cpu().numpy()



'''3. vector operations'''
def dotproduct(V1,V2):
    '''dot product of two vectors'''
    if type(V1.x)==np.ndarray:
        v1=np.concatenate(([V1.x.ravel()],[V1.y.ravel()],[V1.z.ravel()]),axis=0)
        v2=np.concatenate(([V2.x.ravel()],[V2.y.ravel()],[V2.z.ravel()]),axis=0)
        product=np.sum(v1*v2,axis=0)
    elif type(V1.x)==T.Tensor:
        v1=T.concatenate(([V1.x.ravel()],[V1.y.ravel()],[V1.z.ravel()]),axis=0)
        v2=T.concatenate(([V2.x.ravel()],[V2.y.ravel()],[V2.z.ravel()]),axis=0)
        product=T.sum(v1*v2,axis=0)
    return product
def scalarproduct(K,V1):
    '''scalar product between a vector and a constent'''
    V2=copy.copy(V1)
    V2.x=K*V2.x
    V2.y=K*V2.y
    V2.z=K*V2.z    
    return V2
def sumvector(V1,V2):
    ''''''
    V3=copy.copy(V1)
    V3.x=V1.x+V2.x
    V3.y=V1.y+V2.y
    V3.z=V1.z+V2.z
    return V3

