import numpy as np

from coordinate import global_coord
from Vectorpy import vector

class Spherical_Field(vector):
    def __init__(self,center,sizex,sizey,Nx,Ny,Region='far',Type='uv',coord_sys=global_coord,distance=300):
        vector.__init__(self)
        self.coord=coord_sys
        '''
        center: map center in angles;
        sizex: angular size in u or Azimuth;
        sizey: angular size in v or Elevation;
        Nx,Ny are the sampling number;
        if 'Region' is 'near', the field is represented in Cartesian coordinates.

        '''
        Grid_type={'uv':     lambda x,y: (x,y,np.sqrt(1-(x**2+y**2))),
               'EloverAz':lambda x,y: (-np.sin(x)*np.cos(y),np.sin(y),np.cos(x)*np.cos(y))}
        
        self.x,self.y=np.meshgrid(np.linspace(-sizex/2,sizex/2,Nx),np.linspace(-sizey/2,sizey/2,Ny))
        self.x=self.x+center[0]
        self.y=self.y+center[1]
        self.u,self.v,self.z=Grid_type[Type](self.x,self.y) # self.u self.v is the angular in local coordinates.

        if Region.lower()=='far':
            #xyz=np.matmul(self.coord.mat_l_g,xyz=np.append([self.x,self.y],[self.z],axis=0))
            #self.x=xyz[0,:]
            #self.y=xyz[1,:]
            #self.z=xyz[2,:]
            pass
        elif Region.lower()=='near':
            self.x=self.x*distance
            self.y=self.y*distance
            self.z=self.z*distance
            self.x,self.y,self.z=self.coord._toGlobal_coord(self.x,self.y,self.z)