"""
This script offers several functions to convert the define electric objectors,
like frequency, PO anlaysis, and so on,
into text str which can be readable for Ticra Grasp SOFTWARE.
"""

import numpy as np;
import transforms3d;

'''5. write frequency'''
def write_frequency(name,freq_list):
    freq=np.array(freq_list).ravel()
    Str='';
    Str+=name+'  frequency\n(\n';
    Str+='  frequency_list   : sequence(';
    for item in freq:
        Str+=str(item)+' GHz,'
    Str+=')\n)\n\n'
    return Str;
    
        
'''6. write gaussian beam'''
def write_beam(name,freq,coor_sys,taper_angle=11.894,taper=-8,polarisation='linear_x'):
    Str=''
    Str+=name+'  gaussian_beam_pattern\n(\n';
    Str+='  frequency      : ref('+freq+'),\n';
    Str+='  coor_sys       : ref('+coor_sys+'),\n';
    Str+='  taper_angle    : '+str(taper_angle)+',\n';
    Str+='  taper          : '+str(taper)+',\n';
    Str+='  polarisation   : '+polarisation+'\n)\n\n';
    
    return Str;
'''7.  write multi po simulations'''    
def write_multi_po(name,freq,scatterer,outputfile='',method='po_plus_ptd'):
    Str='';
    Str+=name+'.po  '+'po_multi_face_scatterer\n(\n'
    Str+='  frequency      : ref('+freq+'),\n';
    Str+='  scatterer      : ref('+scatterer+'),\n';
    Str+='  method         : '+method+',\n';
    Str+='  ptd_contributions : all,\n';
    if outputfile=='':
        Str+='  file_name      : '+name+'.po.cur\n)\n\n';
    else:
        Str+='  file_name      : '+outputfile+'\n)\n\n';
    return Str;
    
    
'''8. write spherical grid field'''
def write_spherical_grid(name,coor_sys,u_range,v_range,u0,v0,Nu,Nv,near_far='near',near_dist=100,filename=''):
    Str='';
    Str+=name+'.grd  spherical_grid\n(\n';
    Str+='  coor_sys       : ref('+coor_sys+'),\n';
    Str+='  x_range        : struct(start: '+str(u0-u_range/2)+', end: '+str(u0+u_range/2)+', np: '+str(int(Nu))+'),\n';
    Str+='  y_range        : struct(start: '+str(v0-v_range/2)+', end: '+str(v0+v_range/2)+', np: '+str(int(Nv))+'),\n';
    Str+='  near_far       : '+near_far+',\n';
    Str+='  near_dist      : '+str(near_dist)+' m,\n';
    if filename=='':
        Str+='  file_name      : '+name+'.grd\n)\n\n';
    else:
        Str+='  file_name      : '+filename+'\n)\n\n';
    
    return Str;

def write_planar_grid(name,coor_sys,near_dist=0,x_range=[-1,1,11],y_range=[-1,1,11],grid_type='xy',filename=''):
    Str='';
    Str+=name+'.grd  planar_grid\n(\n';
    Str+='  coor_sys        : ref('+coor_sys+'),\n';
    Str+='  near_dist       : '+str(near_dist)+' m,\n';
    Str+='  grid_type       : '+grid_type+',\n';
    Str+='  x_range         : struct(start:'+str(x_range[0])+', end:'+str(x_range[1])+', np:'+str(int(x_range[2]))+', unit: mm),\n';
    Str+='  y_range         : struct(start:'+str(y_range[0])+', end:'+str(y_range[1])+', np:'+str(int(y_range[2]))+'),\n';
    if filename=='':
        Str+='  file_name       : '+name+'.grd,\n';
    else:
        Str+='  file_name       : '+filename+',\n';
        
    Str+=')\n\n';
    return Str;
    




