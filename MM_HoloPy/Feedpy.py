'''
This package provides N input beams, and each beam function can offer scalar and vector modes.
1. Gaussian beam in far field;
2. Gaussian beam near field;
'''

import numpy as np;
from coordinate import global_coord
c=299792458*1000#mm/s _light_speed
'''
Type 1: Gaussian beam;
'''
class Gaussianfeed():
    '''
    feed horn
    '''
    def __init__(self,freq,E_taper=-12,E_angle=13,coord_sys=global_coord,polarization='scalar'):
        '''
        E_taper is the illumination edge taper for the telescope.
        E_angle is the illumination angular size with the pre-defined edge taper. 
        '''
        self.freq=freq*10**9 #GHz
        self.Lambda=c/self.freq #mm
        self.k=np.pi*2/self.Lambda
        
        self.coord=coord_sys
        self.E_taper=E_taper
        self.E_angle=E_angle
        self.polar=polarization

    def field(self,Mirror):
        k=self.k
        M=Mirror.to_coord(self.coord)
        M_n=Mirror.vector_rotation(self.coord)
        r,theta,phi = self.coord._toSpherical(M.x,M.y,M.z)
        if self.polar.lower()=='scalar':
            Theta_max=self.E_taper
            E_taper=self.E_taper
            
            b=(20*np.log10((1+np.cos(Theta_max))/2)-E_taper)/(20*k*(1-np.cos(Theta_max))*np.log10(np.exp(1)))
            w0=np.sqrt(2/k*b)
            
            R=np.sqrt(r**2-b**2+1j*2*b*M.z)
            E=np.exp(-1j*k*R-k*b)/R*(1+np.cos(theta))/2/k/w0*b
            E=E*np.sqrt(8)
                    
            cos_i=np.abs(M.x*M_n.x+M.y*M_n.y+M.z*M_n.z)/r # cosin of reflection angle or incident angle at the reflecting point.

            return E,cos_i