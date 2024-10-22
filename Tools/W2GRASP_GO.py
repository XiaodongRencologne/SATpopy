"""
This script offers several functions to convert the define geometrical elements, reflectors, 
lens, surface, rim and coordinate system into text str which can be readable for Ticra Grasp SOFTWARE.
"""

import numpy as np;
import transforms3d;


# In[2]:


'''1. write coordinate system'''
def write_coord(name,origin,angle,base='mechanical_axis.cs',Type='coor_sys'):
    """
    name, name of the defined coordinate system, coor_reflector1, coor_reflector2....
    origin, the origin of the coordinate in its reference coord_system. 
    
    Normally we define a empty coord system as the first reference, which is called 'M' 
    """
    mat= transforms3d.euler.euler2mat(-angle[0],-angle[1],-angle[2],axes='sxyz')
    x=mat[0,...]
    y=mat[1,...]
    Str=''
    Str+=name+'  '+Type+'\n(\n'
    if name=='mechanical_axis' or name=='mechanical_axis.cs':
        Str+=')\n\n'
    else:
        Str+=' origin     : struct(x: '+str(origin[0])+' mm, y: '+str(origin[1])+' mm, z: '+str(origin[2])+' mm),\n'
        if Type=='coor_sys':            
            Str+=' x_axis     : struct(x: '+str(x[0])+', y: '+str(x[1])+', z: '+str(x[2])+'),\n'
            Str+=' y_axis     : struct(x: '+str(y[0])+', y: '+str(y[1])+', z: '+str(y[2])+'),\n'
            
        if Type=='coor_sys_euler_angles':
            a,b,c=transforms3d.euler.mat2euler(mat,axes='szyz')
            alpha=-a;beta=-b;gamma=-c;            
            Str+=' alpha      : '+str(alpha/np.pi*180)+',\n'
            Str+=' beta       : '+str(beta/np.pi*180)+',\n'
            Str+=' gamma      : '+str(gamma/np.pi*180)+',\n'
        if Type=='coor_sys_grasp_angles':
            a,b,c=transforms3d.euler.mat2euler(mat,axes='szyz')
            phi=-a
            theta=-b
            psi=-c+phi
            Str+=' theta      : '+str(theta/np.pi*180)+',\n'
            Str+=' phi        : '+str(phi/np.pi*180)+',\n'
            Str+=' psi        : '+str(psi/np.pi*180)+',\n'           
        Str+=' base       : ref('+base+')\n'
        Str+=')\n\n' 
    return Str

'''2. write surface'''
def write_Tabulated_surf(name,filename):
    Type='regular_xy_grid'
    Str=''
    Str+=name+'.surf  '+Type+'\n(\n'
    Str+='  file_name : '+filename+',\n'
    Str+='  xy_unit   : mm,\n'
    Str+='  z_unit    : mm\n)\n'
    return Str

def write_Tabulated_Rota_symm_surf(name,filename):
    Type='rotationally_symmetric'
    Str=''
    Str+=name+'.surf  '+Type+'\n(\n'
    Str+='  file_name : '+filename+',\n'
    Str+='  rho_unit   : mm,\n'
    Str+='  z_unit    : mm\n)\n'
    return Str

def write_plane_surf(name,point):
    Str=''
    Str+=name+'  plane\n(\n'
    Str+='point       : struct(x: '+str(point[0])+' mm, y: '+str(point[1])+' mm, z: '+str(point[2])+' mm)\n)\n\n'
    return Str

def write_conic_surf(name,r1,r2,theta_i,theta_n):
    '''r1 r2 is the distance between F1 or F2 and reflection point, with unit;
       theta_i reflection angle
       theta_n normal vector angle from axis z to n, positive around y-axis.
    '''
    Str = ''
    Str += name + '  conical_surface\n(\n'
    Str +='  r1      : '+ r1 +',\n'
    Str +='  r2      : '+ r2 +',\n'
    Str +='  theta_i : '+ str(theta_i) +',\n'
    Str +='  theta_n : '+ str(theta_n) +'\n)\n'



'''3. write rim'''
def write_rim(name,centre,side_lengths,Type='rectangular_rim'):
    """
    Type: rectangular_rim
          elliptical_rim
    """
    Str=''
    Str+=name+'.rim  '+Type+'\n(\n'
    Str+='  centre       : struct(x: '+str(centre[0])+' mm, y: '+str(centre[1])+' mm),\n'
    if Type == 'rectangular_rim':
        Str+='  side_lengths : struct(x: '+str(side_lengths[0])+' mm, y: '+str(side_lengths[1])+' mm)\n)\n'
    elif Type == 'elliptical_rim':
        Str+='  half_axis    : struct(x: '+str(side_lengths[0])+' mm, y: '+str(side_lengths[1])+' mm)\n)\n'
    else:
        raise ValueError("input rim type is wrong!!!")
    return Str

'''4. write scatter'''

def write_refl_multiple(name,coord_sys,surfs,rims,distortions=''):
    Str=''
    Str+=name+'.refl  '+'panels_individually_defined\n(\n'
    Str+='  coor_sys      : ref('+''+coord_sys+'),\n'
    Str+='  panels        : sequence\n    ('
    n=0
    for item in rims:
        if distortions=='':
            Str+='    struct(surface: ref('+surfs+'), rim: ref('+item+'),distortion: ref('+distortions+'),hinge_coor_sys: ref(), hinge_rotation: 0.0),\n';
        else:
            Str+='    struct(surface: ref('+surfs+'), rim: ref('+item+'),distortion: ref('+distortions[n]+'),hinge_coor_sys: ref(), hinge_rotation: 0.0),\n'
            n+=1
    Str+='    )\n'
    Str+=')\n'
    
    return Str


'''5. write a simple lens'''

class coating_surface():
    def __init__(self, thickness, 
                 refractive_index,
                 loss_tangent=0):
        self.t0 = thickness
        self.r_index = str(refractive_index)
        self.loss_t = str(loss_tangent)

        self.Str = 'struct(thickness:' + self.t0+',\n'
        self.Str += 'refractive_index:' +self.r_index+',\n'
        self.Str += 'loss_tangent:' +self.loss_t+')'

def write_simple_lens(name,coord_sys,diameter,
                      refractive_index, loss_tangent=0,
                      r1=0, r2=0,
                      bs1=0,bs2=0,
                      thickness=0,
                      surf1_file='',surf2_file='',lengthUnit_file='mm',
                      coating_surf1='',coating_surf2=''):
    Str = ''
    Str += name + 'simple_lens\n(\n'
    Str += '  coor_sys         :ref('+''+coord_sys+'),\n'
    Str += '  diameter         :' + diameter+',\n'
    Str += '  refractive_index :' + str(refractive_index) +',\n'
    Str += '  loss_tangent     :' + str(loss_tangent) +',\n'
    if r1 != 0:
        Str += '     r1           :'+ r1 +',\n'
    if r2 != 0:
        Str += '     r2           :'+ r2 +',\n'
    Str += '  bs1              :' + str(bs1) +',\n'
    Str += '  bs2              :' + str(bs2) +',\n'
    Str += '  thickness        :' + str(thickness) +',\n'
    if surf1_file != '':
        Str += '  surface1_file    :' + str(surf1_file) +',\n'
    if surf2_file != '':
        Str += '  surface2_file       :' + str(surf2_file) +',\n'
    Str +=     '  length_unit_in_files :' + lengthUnit_file +',\n'
    if coating_surf1 != '':
        Str += '  coating_surface1     :' + coating_surf1 +',\n'
    if coating_surf2 != '':
        Str += '  coating_surface2     :' + coating_surf2 +',\n' 

    Str += ')\n'
    return Str
    



'''6.  write aperture in screen'''

def write_aperture_scan(name,coord_sys,rim,infinity_shadow='on'):
    Str = ''
    Str += name + 'aperture_in_screen\n(\n'
    Str += '  coor_sys      : ref(' + coord_sys+'),\n'
    Str += '  rim           : ref(' + rim +'),\n'
    Str += '  infinity_shadow: '+infinity_shadow +'\n'
    Str += ')\n'
    return Str