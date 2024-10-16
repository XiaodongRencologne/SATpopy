# %%
import numpy as np
import matplotlib.pyplot as plt
import copy

import pyvista as pv
#pv.set_jupyter_backend('trame')
#pv.set_jupyter_backend(None) 
#pv.set_jupyter_backend('server')
#pv.set_jupyter_backend('static')

from coordinate import global_coord
from rim import elliptical_rim,rect_rim
from surface import PolySurf
from Vectorpy import vector


p=pv.Plotter()

# %%
class Geo():
    '''
    Basic optical element.
    '''
    def __init__(self,surf,rim,coord_sys=global_coord):
        '''basic parameters'''
        self.surf=surf
        self.rim=rim
        self.coord=coord_sys
        self.View_3D=None
        '''coordinates of the object'''
        self.points=vector()
        self.r_points=vector()
        self.g_points=vector()
        self.normal_v=vector()
        

    def view(self,widget):
        x,y,w=self.rim.sampling(10,10)
        del(w)
        z=self.surf.surface(x,y)
        x,y,z=self.coord._toGlobal_coord(x,y,z)
        points = np.c_[x.reshape(-1), y.reshape(-1), z.reshape(-1)]
        del(x,y,z)
        if self.View_3D!=None:
            widget.remove_actor(self.View_3D)
        else:
            cloud = pv.PolyData(points)
            self.View_3D = cloud.delaunay_2d()
        del(points)
        widget.add_mesh(self.View_3D,show_edges=True)
        
    def sampling(self,Nx,Ny):
        # sampling
        self.points.x,self.points.y,self.weight=self.rim.sampling(Nx,Ny)
        # calculate surface in z.
        self.points.z=self.surf.surface(self.points.x,self.points.y)
        # to global coordinates
        self.g_points.x,self.g_points.y,self.g_points.z=self.coord._toGlobal_coord(self.g_points.x,self.g_points.y,self.g_points.z)
        # surface normal vectors.
        self.normal_v.x,self.normal_v.y,self.normal_v.z,self.N=self.surf.normal_vector(self.points.x,self.points.y)
        xyz=np.append([self.normal_v.x,self.normal_v.y],[self.normal_v.z],axis=0)
        # convert normal vector in global coordinate system
        self.normal_v.x,self.normal_v.y,self.normal_v.z=np.matmul(self.coord.mat_l_g,xyz)
        del(xyz)

    def to_coord(self,tar_coord):
        points_new=vector()
        points_new.x,points_new.y,points_new.z=tar_coord.Global_to_local(self.g_points.x,self.g_points.y,self.g_points.z)
        return points_new
    
    def vector_rotation(self,tar_coord):
        Nvector=copy.copy(self.normal_v)
        xyz=np.append([self.normal_v.x,self.normal_v.y],[self.normal_v.z],axis=0)
        xyz=np.matmul(tar_coord.mat_g_l,xyz)
        Nvector.x=xyz[0,:]
        Nvector.y=xyz[1,:]
        Nvector.z=xyz[2,:]
        return Nvector


# %%
class Reflector(Geo):
    '''
    Reflector, perfect conductor
    '''
    def __init__(self,surf,rim,coord_sys=global_coord,loss=None):
        Geo.__init__(self,surf,rim,coord_sys)