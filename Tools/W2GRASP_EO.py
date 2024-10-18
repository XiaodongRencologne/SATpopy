"""
This script offers several functions to convert the define electric objectors,
like frequency, PO anlaysis, and so on,
into text str which can be readable for Ticra Grasp SOFTWARE.
"""

import numpy as np;
import transforms3d;

'''1. write frequency list'''
def write_frequency_list(name,freq_list):
    freq=np.array(freq_list).ravel()
    Str=''
    Str+=name+'  frequency\n(\n'
    Str+='  frequency_list   : sequence('
    for item in freq:
        Str+=str(item)+' GHz,'
    Str+=')\n)\n\n'
    return Str
    
        
'''2. write  Feed: pattern: gaussian beam'''
def write_Gauss_beam(name,freq,
                     coor_sys,
                     taper_angle=11.894,taper=-8,
                     polarisation='linear_x'):
    Str=''
    Str+=name+'  gaussian_beam_pattern\n(\n'
    Str+='  frequency      : ref('+freq+'),\n'
    Str+='  coor_sys       : ref('+coor_sys+'),\n'
    Str+='  taper_angle    : '+str(taper_angle)+',\n'
    Str+='  taper          : '+str(taper)+',\n'
    Str+='  polarisation   : '+polarisation+'\n)\n\n'
    return Str

'''3. write  Feed: pattern: gaussian beam in near field or using Qausi optics beam'''
def write_Gauss_beam_near(name,freq,
                          coor_sys,
                          beam_radius, phase_front_radius,
                          polarisation='linear_x',
                          factor = [0,0], frequency_index_for_plot = 1):
    """
    beam_radius: w = w_0 sqrt(1+z0^2/b^2)
    phase front radius, Rx = z0(1+b^2/z0^2)
    b = w0^2k/2
    factor is a multiplier [Amp(dB),phase(deg)]
    """
    Str=''
    Str+=name+'  gaussian_beam_pattern\n(\n'
    Str+='  frequency      : ref('+freq+'),\n'
    Str+='  coor_sys       : ref('+coor_sys+'),\n'
    Str+='  beam_radius    : '+ beam_radius+',\n'
    Str+='  phase_front_radius : '+phase_front_radius+',\n'
    Str+='  polarisation   : '+polarisation+',\n'
    Str+='  factor         :struct(db:'+str(factor[0])+',deg:'+str(factor[1])+'),\n'
    Str+='  frequency_index_for_plot :'+str(frequency_index_for_plot)+'\n)\n'
    return Str

### PO analysis
'''4. write PO single Face scatterer'''
def write_single_po(name,freq,
                    scatterer,
                    method=0,
                    po_points = [0,0],
                    ptd_points= [[-1,0]],
                    factor= [0,0],
                    spill_over='off',
                    coord_sys=None,
                    outputfile = ''):
    Method=['po_plus_ptd',
            'po',
            'ptd']
    Str=''
    Str+= name +'.po  '+'   po_single_face_scatterer\n(\n'
    Str+='  frequency      :ref('+freq+'),\n'
    Str+='  scatterer      :ref('+scatterer+'),\n'
    Str+='  method         :'+Method[method]+',\n'
    Str+='  po_points      :'+'struct(po1:'+str(po_points[0])+', po2:' + str(po_points[1])+'),\n'
    Str+='  ptd_points     :'+'sequence(\n'
    for item in ptd_points:
        Str+='                              struct(edge:'+ str(item[0])+',\n'
        Str+='                                     ptd:'+ str(item[1])+'),\n'
    Str+='                              ),\n'
    Str+='  factor         :struct(db:'+str(factor[0])+', deg:'+str(factor[1])+'),\n'
    Str+='  spill_over     :'+spill_over+',\n'
    if coord_sys !=None:
        Str+='  coor_sys       :ref('+coord_sys+'),\n'
    if outputfile =='':
        Str+='  file_name  :'+name+'.po.cur\n)\n\n'
    else:
        Str+='  file_name      : '+outputfile+'\n)\n\n'
    return Str

'''4.  write multi po simulations'''    
def write_multi_po(name,freq,
                   scatterer,
                   outputfile='',
                   method=0):
    Method=['po_plus_ptd',
            'po',
            'ptd']
    Str=''
    Str+=name+'.po  '+'po_multi_face_scatterer\n(\n'
    Str+='  frequency      : ref('+freq+'),\n'
    Str+='  scatterer      : ref('+scatterer+'),\n'
    Str+='  method         : '+Method[method]+',\n'
    Str+='  ptd_contributions : all,\n'
    if outputfile=='':
        Str+='  file_name      : '+name+'.po.cur\n)\n\n'
    else:
        Str+='  file_name      : '+outputfile+'\n)\n\n'
    return Str
    
'''5. write PO lens'''

def write_lens_po(name, freq, 
                  lens, 
                  get_field='lens_in_screen',
                  method='go_plus_po', waist_radius=0,
                  po_points={'face1':[0,0],'face2': [0,0]}, factor=[0,0],
                  spill_over='off',coord_sys='', 
                  current_file_face1='', 
                  current_file_face2='',
                  gbc_file=''):
    '''
    get field: lens_in_screen, lens_in_free_space, face1, face2
    
    '''
    Str=''
    Str+= name + '.po  '+ 'po_lens\n(\n'
    Str+= '  frequence     :ref('+freq+'),\n'
    Str+= '  lens          :ref('+lens+'),\n'
    Str+= '  get_field     :'+get_field+',\n'
    Str+= '  method        :'+method+',\n'
    if method == 'gauss_laguerre':
        Str+= '  method        :'+method+',\n'
        if waist_radius==0:
            Str+= '  waist_radius  :0,\n'
        else:
            Str+= '  waist_radius  :squence(\n'
            for item in waist_radius:
                Str+= '  waist_radius  :squence(\n'
                Str+= '                        '+str(item)+',\n'
            Str+='                          ),\n'
    Str+= '  po_points     :struct(face1_po1:'+str(po_points['face1'][0])
    Str+= ', face1_po2:'+str(po_points['face1'][1])+',\n'
    Str+= '                           face2_po1:'+str(po_points['face2'][0])
    Str+= ', face2_po2:'+str(po_points['face2'][1])+'),\n'
    Str+= '  factor        :struct(db:'+str(factor[0])+',deg:'+str(factor[1])+'),\n'
    Str+= 'spill_over      :'+spill_over+',\n'
    if coord_sys !='':
        Str+= 'coor_sys        :ref('+coord_sys+'),\n'
    if current_file_face1 =='':
        Str+= 'current_file_face1 :'+name+'face1.po.cur,\n'
    else:
        Str+= 'current_file_face1 :'+current_file_face1+'.po.cur,\n'
    if current_file_face2 =='':
        Str+= 'current_file_face2 :'+name+'face2.po.cur,\n'                        
    else:
        Str+= 'current_file_face2 :'+current_file_face2+'.po.cur,\n'

    Str+=')\n\n'

''' write Aperture in screen PO'''

#MoM

''' write MoM'''
def write_MoM(name,):
    pass

'''8. write spherical grid field'''
def write_spherical_grid(name,coor_sys,u_range,v_range,u0,v0,Nu,Nv,near_far='near',near_dist=100,filename=''):
    Str=''
    Str+=name+'.grd  spherical_grid\n(\n'
    Str+='  coor_sys       : ref('+coor_sys+'),\n'
    Str+='  x_range        : struct(start: '+str(u0-u_range/2)+', end: '+str(u0+u_range/2)+', np: '+str(int(Nu))+'),\n'
    Str+='  y_range        : struct(start: '+str(v0-v_range/2)+', end: '+str(v0+v_range/2)+', np: '+str(int(Nv))+'),\n'
    Str+='  near_far       : '+near_far+',\n'
    Str+='  near_dist      : '+str(near_dist)+' m,\n'
    if filename=='':
        Str+='  file_name      : '+name+'.grd\n)\n\n'
    else:
        Str+='  file_name      : '+filename+'\n)\n\n'
    
    return Str

def write_planar_grid(name,coor_sys,near_dist=0,x_range=[-1,1,11],y_range=[-1,1,11],grid_type='xy',filename=''):
    Str=''
    Str+=name+'.grd  planar_grid\n(\n'
    Str+='  coor_sys        : ref('+coor_sys+'),\n'
    Str+='  near_dist       : '+str(near_dist)+' m,\n'
    Str+='  grid_type       : '+grid_type+',\n'
    Str+='  x_range         : struct(start:'+str(x_range[0])+', end:'+str(x_range[1])+', np:'+str(int(x_range[2]))+', unit: mm),\n'
    Str+='  y_range         : struct(start:'+str(y_range[0])+', end:'+str(y_range[1])+', np:'+str(int(y_range[2]))+'),\n'
    if filename=='':
        Str+='  file_name       : '+name+'.grd,\n'
    else:
        Str+='  file_name       : '+filename+',\n'
        
    Str+=')\n\n'
    return Str
    
