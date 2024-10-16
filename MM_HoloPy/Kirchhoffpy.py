#%%
import numpy as np;
import torch as T;
import copy;
import time;
from multiprocessing import Pool
# %%

'''1. Define Kirchhoff-Fresnel Integration Solver.'''
#def KirchhoffSolver(Source, Tar_element, Type='real', Field_type='far'):
def kirchhoff(m1,m2,Field_m1,cos_i,k,Keepmatrix=False,Device='cpu',**keywards):
    '''
    Formulatar of the Kirchhoff-Fresnel Integration method
    *******Field=Sum{field()*exp(-j*k*r)/R*[cos_i+cos_r] J*dx*dy}*************
    **************************************************************************
    m1: reflector of source field;
    m2: target reflector;
    Field_in: Field on 'm1';
    cos_i: cosin of reflection angle on m1, in the following document 
    cos_r is cosin of the output ray from the point on 'm1' to target point on 'm2' ;
    Keepmatrix means if the intermediate computing results are saved as a 2D matrix used for speeding up
    calculations.'
    Device: 'cpu' means the computation is done by cpu. 'cuda' can load the data into
    GPU RAM and computed by GPU for acceleration.
    '''
    m1.g_points.to_Tensor(Device)
    m2.g_points.to_Tensor(Device)
    cos_i=T.tensor(cos_i,dtype=T.float64).to(T.device(Device))
    cos_i = cos_i.to(Device)
    Field_m1=T.tensor(Field_m1,dtype=T.complex128).to(Device)
    if Device=='cuda':
        print('The computation is speed up by GPU units!')            
    # Define output field:
    #Field_m2=T.tensor()
    def calcu(x,y,z,F_m1):
        x=m2.g_points.x.reshape(-1,1)-m1.g_points.x.reshape(1,-1)
        y=m2.g_points.y.reshape(-1,1)-m1.g_points.y.reshape(1,-1)
        z=m2.g_points.z.reshape(-1,1)-m1.g_points.z.reshape(1,-1)
        #cosin terms
        if cos_i.size==1:
            cos=1
            x=T.sqrt(x**2+y**2+z**2)
        else:
            cos=T.abs(x*m1.normal_v.x.reshape(1,-1)
                        +y*m1.normal_v.y.reshape(1,-1)
                        +z*m1.normal_v.z.reshape(1,-1))
            x=T.sqrt(x**2+y**2+z**2)
            cos=(cos/x+np.abs(cos_i))/2
        del(y,z)
        # field calculation    
        Amp=1/x*m1.N*m1.weight/2/T.pi*k*cos #
        Mat0=Amp*np.exp(-1j*k*x)
        F_m2=T.sum(Mat0*F_m1.reshape(1,-1),axis=-1)
        return Mat0, F_m2
    