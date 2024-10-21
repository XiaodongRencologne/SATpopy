import numpy as np
import W2GRASP_GO, W2GRASP_EO
from W2GRASP_GO import write_coord,write_simple_lens


'''
Define Geometrical objectors, coordinate system, simple lens.
'''
class coor_sys():
    def __init__(self,origin, angle, ref_coor_name = 'mechanical_axis.cs', name='coor_sys'):
        self.name = name
        self.ref_name = ref_coor_name
        self.origin = origin
        self.angle = angle

        self.Str = write_coord(self.name,
                               self.origin,
                               self.angle,
                               base = self.ref_name)
        self.List = {'name': 'coor_feed_ref',
                     'origin':[0,0,0],
                     'angle': [0,0,0],
                     'base':'mechanical_axis.cs',
                     'Type':'coor_sys'}
        

class simple_lens():
    def __init__(self, coor_sys, 
                 diameter, 
                 refractive_index, loss_tangent = 0, 
                 r1 = None, r2 = None, bs1 = 0, bs2 = 0,
                 thickness = 0,
                 surf_f1 = '',
                 surf_f2 = '',
                 lengthUnit = 'mm',
                 coating_surf1 = '', coating_surf2 = '',
                 name='simple_lens'):
        self.name = name
        self.coor_name = coor_sys
        self.diameter = diameter
        self.refractive_index = refractive_index
        self.loss_tangent = loss_tangent
        self.r1 =r1
        self.r2 = r2
        self.bs1 = bs1
        self.bs2 = bs2
        self.thickness =thickness
        self.surf_f1 = surf_f1
        self.surf_f2 = surf_f2
        self.lengthUnit = lengthUnit
        self.coating_surf1 = coating_surf1
        self.coating_surf2 = coating_surf2

        self.Str = write_simple_lens(self.name,
                                     self.coor_name,
                                     str(self.diameter)+'mm',
                                     self.refractive_index,
                                     self.loss_tangent,
                                     r1 = self.r1, r2 = self.r2, 
                                     bs1 = self.bs1, bs2 = self.bs2,
                                     thickness = self.thickness,
                                     surf1_file = self.surf_f1,surf2_file = self.surf_f2,
                                     lengthUnit_file = self.lengthUnit,
                                     coating_surf1 = self.coating_surf1, coating_surf2 = self.coating_surf2)