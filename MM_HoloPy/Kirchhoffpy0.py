#%%
import numpy as np;
import torch as T;
from Vectorpy import Fvector
import copy;
import time;
# %%

'''1. Define Kirchhoff-Fresnel Integration Solver.'''
#def KirchhoffSolver(Source, Tar_element, Type='real', Field_type='far'):
def kirchhoff(m1,m2,Field_m1,cos_i,k,Keepmatrix=False,**keywards):
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
    Status_device='cpu'
    #load the data into Tensor
    try:
        m1.g_points.to_Tensor('cuda')
        m2.g_points.to_Tensor('cuda')
        cos_i=T.tensor(cos_i,dtype=T.float64).to(T.device('cuda'))
        Field_m1.to_Tensor('cuda')
        Status_device='cuda'
        print('The computation is speed up by GPU units!')
    except:
        m1.g_points.to_Tensor('cpu')
        m2.g_points.to_Tensor('cpu')
        cos_i=T.tensor(cos_i,dtype=T.float64).to(T.device('cpu'))
        Field_m1.to_Tensor('cpu')
        Status_device='cpu'

    # Define output field:
    Field_m2=Fvector(m2.g_points.x.size)
    cos_i=cos_i.reshape(1,-1)
    




    if Device.lower()=='cpu':
        m1.g_points.to_numpy()
        m2.g_points.to_numpy()
        if Keepmatrix:
            Mat=np.array([],dtype=np.complex128)
            for i in range(m2.g_points.x.size):
                x=m2.g_points.x[i]-m1.g_points.x.reshape(1,-1)
                y=m2.g_points.y[i]-m1.g_points.y.reshape(1,-1)
                z=m2.g_points.z[i]-m1.g_points.z.reshape(1,-1)
                r=np.sqrt(x**2+y**2+z**2)
                
                #cosin terms
                cos_r=np.abs(x*m1.normal_v.x.reshape(1,-1)
                            +y*m1.normal_v.y.reshape(1,-1)
                            +z*m1.normal_v.z.reshape(1,-1))/r
                cos=(np.abs(cos_r)+np.abs(cos_i))/2
                if cos_i.size==1:
                    cos=1
                # field calculation    
                Amp=1/r*m1.N*m1.weight/2/np.pi*np.abs(k)*cos
                phase=-k*r

                Mat0=Amp*np.exp(1j*phase)
                Field[i]=Mat0*Field_in.reshape(1,-1)
                if i==int(m2.g_points.x.size/2):
                    COS_R=cos_r
                Mat=np.concatenate((Mat,Mat0),axis=0)
        else:
            for i in range(m2.g_points.x.size):
                x=m2.g_points.x[i]-m1.g_points.x.reshape(1,-1)
                y=m2.g_points.y[i]-m1.g_points.y.reshape(1,-1)
                z=m2.g_points.z[i]-m1.g_points.z.reshape(1,-1)
                r=np.sqrt(x**2+y**2+z**2)
                
                #cosin terms
                cos_r=np.abs(x*m1.normal_v.x.reshape(1,-1)
                             +y*m1.normal_v.y.reshape(1,-1)
                             +z*m1.normal_v.z.reshape(1,-1))/r
                cos=(np.abs(cos_r)+np.abs(cos_i))/2
                if cos_i.size==1:
                    cos=1
                # field calculation    
                Amp=1/r*m1.N*m1.weight/2/np.pi*np.abs(k)*cos
                phase=-k*r

                Mat=Amp*np.exp(1j*phase)
                Field[i]=Mat*Field_in.reshape(1,-1)
                if i==int(m2.g_points.x.size/2):
                    COS_R=cos_r
    # calculation in GPU.
    elif Device.lower()=='cuda':
        m1.g_points.to_Tensor('cuda')
        m2.g_points.to_Tensor('cuda')
        cos_i=T.tensor(cos_i,dtype=T.float64).to(T.device('cuda'))
        Field_in=T.tensor(Field_in,dtype=T.cdouble).to(T.device('cuda'))
        if Keepmatrix:
            Mat=T.tensor(np.array([]),dtype=T.cdouble)
            for i in range(m2.g_points.x.size):
                x=m2.g_points.x[i]-m1.g_points.x.reshape(1,-1)
                y=m2.g_points.y[i]-m1.g_points.y.reshape(1,-1)
                z=m2.g_points.z[i]-m1.g_points.z.reshape(1,-1)
                r=T.sqrt(x**2+y**2+z**2)
                
                #cosin terms
                cos_r=T.abs(x*m1.normal_v.x.reshape(1,-1)
                             +y*m1.normal_v.y.reshape(1,-1)
                             +z*m1.normal_v.z.reshape(1,-1))/r
                cos=(T.abs(cos_r)+T.abs(cos_i))/2
                if cos_i.size==1:
                    cos=1
                # field calculation    
                Amp=1/r*m1.N*m1.weight/2/T.pi*T.abs(k)*cos
                phase=-k*r

                Mat0=Amp*(T.cos(phase)+1j*T.sin(phase))
                Field[i]=Mat0*Field_in.reshape(1,-1)
                Mat=T.cat((Mat,Mat0),axis=0)
                if i==int(m2.g_points.x.size/2):
                    COS_R=cos_r
            Mat=Mat.reshape(m2.g_points.x.size,-1)
        else:
            for i in range(m2.g_points.x.size):
                x=m2.g_points.x[i]-m1.g_points.x.reshape(1,-1)
                y=m2.g_points.y[i]-m1.g_points.y.reshape(1,-1)
                z=m2.g_points.z[i]-m1.g_points.z.reshape(1,-1)
                r=T.sqrt(x**2+y**2+z**2)
                
                #cosin terms
                cos_r=T.abs(x*m1.normal_v.x.reshape(1,-1)
                             +y*m1.normal_v.y.reshape(1,-1)
                             +z*m1.normal_v.z.reshape(1,-1))/r
                cos=(T.abs(cos_r)+T.abs(cos_i))/2
                if cos_i.size==1:
                    cos=1
                # field calculation    
                Amp=1/r*m1.N*m1.weight/2/T.pi*T.abs(k)*cos
                phase=-k*r

                Mat=Amp*(T.cos(phase)+1j*T.sin(phase))
                Field[i]=Mat0*Field_in.reshape(1,-1)
                if i==int(m2.g_points.x.size/2):
                    COS_R=cos_r
    return Mat, Field, COS_R