# %%
import os
import re

import numpy as np
import scipy
from scipy import interpolate

# %%
def W_surf(x0, x1, y0, y1, Nx, Ny, z_xy, filename):
    '''create a spline surface file!!!'''
    header='x,  '+str(x0)+', ' + str(x1)+'\n'\
           + 'Nx, '+str(Nx)+'\n'\
           + 'y,  '+str(y0)+', ' + str(y1)+'\n'\
           + 'Ny, '+str(Ny)
    
    if z_xy.size==Nx*Ny:
        np.savetxt(filename+'.surf', z_xy.ravel(), header= header, delimiter=',')
        return True
    else:
        print('points number does not agree with each other!')
        return False

def R_surf(surf_file):
    with open(surf_file,'r') as f:
        line=f.readline()
        string=re.split(',| |:',line.replace(" ", ""))
        x0=float(string[1])
        x1=float(string[2])

        line=f.readline()
        string=re.split(',| |:',line.replace(" ", ""))
        Nx=int(string[1])

        line=f.readline()
        string=re.split(',| |:',line.replace(" ", ""))
        y0=float(string[1])
        y1=float(string[2])

        line=f.readline()
        string=re.split(',| |:',line.replace(" ", ""))
        Ny=int(string[1]) 
        f.close()
    
    z=np.genfromtxt(surf_file,delimiter=',',skip_header=4)

    return x0, x1, Nx, y0, y1, Ny, z

### write zemax symmetrical lens surface
def zemax2RSF(Np,Kspace,Ktip,lens_para,outputfolder='',sign = 1):
    '''
    Rotationally symmetric surface.
    It is one-demensional surface which is a function of radial Rho. 
    rho =sqrt(x^2+y^2)

    Lens_para = {'R': 500,
                 'K': -2.1,
                 'type': 'EvenAsphere',
                 'co': [1,2,3],
                 'D' : 200,
                 'name':'lens_face1'}
    '''
    with open(outputfolder+lens_para['name']+'.rsf','w') as f:
        f.writelines(lens_para['name']+'\n')
        f.writelines(str(Np)+' '+str(Kspace)+' '+str(Ktip)+'\n')
    D = lens_para['r']*2
    if lens_para['type'] == 'EvenAsphere':
        surf_fuc = EvenAsphere(lens_para['R'],lens_para['K'],
                               lens_para['co'])
    if Kspace == 1:
        rho = np.linspace(0,D/2,Np)
        z = sign * surf_fuc(rho)
        data = np.append(rho,z).reshape(2,-1).T
        with open(outputfolder+lens_para['name'] + '.rsf','a') as f:
            #f.writelines(str(rho.min())+' '+str(rho.max()) +'\n')
            np.savetxt(f,data,delimiter=' ',fmt = '%10.8f')
            """
            for n in range(Np):
                f.writelines(str(rho[n]) + ' ' +str(z[n])+'\n')
            """
    return rho, z
def R_lens_surf(surf_file):
    data = np.genfromtxt(surf_file, delimiter=' ', skip_header =2)
    rho = data[:,0]
    z = data[:,1]
    return rho, z
# %%
class PolySurf():
    '''
    Define a surface described by 2-D polynomials.
    
    2-D surface is defined by the polynomials: 
    
      p(x,y) = SUM_ij c_ij * (x/R)^i * (y/R)^j, 
    
    *coefficients are represented by C_ij matrix
    *R is the normalization factor.

    '''
    def __init__(self, coefficients, R=1):
        self.R=R
        if isinstance(coefficients,np.ndarray) or isinstance(coefficients, list):
            self.coefficients=np.array(coefficients)
        elif isinstance(coefficients, str):
            if coefficients.split('.')[-1].lower()=='surfc':
                self.coefficients=np.genfromtxt(coefficients,delimiter=',')
            else:
                print('Please give correct surface coefficients files!')
        else:
            print('The input coefficient list or numpy.ndarry is incorrect!')
        

    def surface(self,x,y):
        z=np.polynomial.polynomial.polyval2d(x/self.R,y/self.R,self.coefficients)
        return z

    def normal_vector(self,x,y,):
        '''normal vector of the surface'''
        nz=-np.ones(x.shape)
        a=np.arange(self.coefficients.shape[0])
        c=self.coefficients*a.reshape(-1,1)
        nx=np.polynomial.polynomial.polyval2d(x/self.R,y/self.R,c[1:,:])/self.R

        a=np.arange(self.coefficients.shape[1])
        c=self.coefficients*a
        ny=np.polynomial.polynomial.polyval2d(x/self.R,y/self.R,c[:,1:])/self.R
        N=np.sqrt(nx**2+ny**2+nz**2)

        nx=nx/N
        ny=ny/N
        nz=nz/N
        '''J: Jacobian Matrix determinant. J=N'''
        return nx,ny,nz,N.ravel()

# %%
class Splines_Surf():
    '''
    Define a surface described by interpolating the input surface data!
    '''
    def __init__(self,surf_file):
        x0,x1,Nx,y0,y1,Ny,z=R_surf(surf_file)
        x=np.linspace(x0,x1,Nx)
        y=np.linspace(y0,y1,Ny)
        
        self._func2d=interpolate.RectBivariateSpline(x,y,z.reshape(Ny,Nx).T)

    def surface(self,x,y):
        z=self._func2d(x,y,grid=False)
        return z
    
    def normal_vector(self,x,y):
        nz = -np.ones(x.shape)
        pass


class Symetrical_surf():
    '''
    Define a rotaional symmetrical surfaces
    '''
    def __init__(self,surf_file):
        rho, z = R_lens_surf(surf_file)
        #self._func1d = interpolate.interp1d(rho, z,kind='cubic')
        self._func1d =interpolate.CubicSpline(rho, z)
    
    def surface(self,x,y):
        return self._func1d(np.sqrt(x**2+y**2))
    
    def normal_vector(self,x,y):
        rho = np.sqrt(x**2+y**2)
        nz = -np.ones(x.shape)
        nx = self._func1d.derivative(1)(rho)* x/rho
        ny = self._func1d.derivative(1)(rho)* y/rho
        N=np.sqrt(nx**2+ny**2+nz**2)

        nx = nx/N
        ny = ny/N
        nz = nz/N

        return nx,ny,nz, N
