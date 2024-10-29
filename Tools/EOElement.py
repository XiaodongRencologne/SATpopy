#%%
import numpy as np
import W2GRASP_GO, W2GRASP_EO
from W2GRASP_EO import write_frequency_list, write_Gauss_beam, write_Gauss_beam_near,write_Gauss_Ellip_Beam
from W2GRASP_EO import write_lens_po, write_Aperture_PO
from W2GRASP_EO import write_spherical_grid, write_spherical_cut
#%%

'''
Define frequency, PO analysis, and MoM analysis
'''
class frequencyList():
    def __init__(self,freq_List,name='freq_range'):
        self.name = name
        self.freq = freq_List
        self.Str = write_frequency_list(self.name, freq_List)


class GaussBeam():
    def __init__(self,freqList,
                 coor_sys,
                 T_angle,Taper,
                 polarisation='linear_x',
                 name='GaussFeed1'):
        self.name = name
        self.freqList = freqList
        self.coor_sys =coor_sys
        self.polarisation =polarisation
        self.T_angle = T_angle
        self.Taper = Taper

        self.Str = write_Gauss_beam(self.name,
                                    self.freqList.name,
                                    self.coor_sys.name,
                                    polarisation = self.polarisation,
                                    taper_angle = self.T_angle,
                                    taper = self.Taper)
        
class GaussBeam_Near():
    def __init__(self,
                 freqList,
                 coor_sys,
                 beam_radius,phase_front_radius,
                 polarisation='linear_x',
                 factor = [0,0],
                 freq_index = 1,
                 name='GaussFeed1'):
        self.coor_sys =coor_sys
        self.freqList = freqList
        self.name = name
        self.b_radius = beam_radius
        self.phase_radius = phase_front_radius
        self.polarisation = polarisation
        self.factor =factor
        self.freq_index =freq_index
        self.Str = write_Gauss_beam_near(self.name,self.freqList.name,
                                         self.coor_sys.name,
                                         self.b_radius, self.phase_radius,
                                         polarisation = self.polarisation,
                                         factor = self.factor,
                                         frequency_index_for_plot = self.freq_index)

class Elliptical_Beam():
    def __init__(self,
                 freqList,
                 coor_sys,
                 Taper, T_angle, 
                 polarisation='linear', 
                 polarisation_angle=0,
                 far_forced = 'off',
                 factor = [0,0],
                 frequency_index_for_plot = 1,
                 name = 'Gaussian_Elliptical_Beam'):
        self.name = name
        self.coor_sys =coor_sys
        self.freqList = freqList
        self.Taper = Taper
        self.T_angle = T_angle
        self.polarisation = polarisation
        self.polar_angle =polarisation_angle
        self.far_forced = far_forced
        self.factor = factor
        self.freq_index = frequency_index_for_plot

        self.Str = write_Gauss_Ellip_Beam(self.name,
                                          self.freqList.name,
                                          self.coor_sys.name,
                                          self.Taper,self.T_angle,
                                          polarisation = self.polarisation,
                                          polarisation_angle = self.polar_angle,
                                          far_forced = self.far_forced,
                                          factor =self.factor,
                                          frequency_index_plot=self.freq_index)



class single_PO():
    def __init__(self):
        pass

class Multi_PO():
    def __init__(self):
        pass

class lens_PO():
    def __init__(self,
                 freq,
                 lens,
                 get_field='lens_in_screen',
                 method='go_plus_po', waist_radius=0,
                 po_points={'face1':[0,0],'face2': [0,0]}, factor=[0,0],
                 spill_over='off',coor_sys='', 
                 current_file_face1='', 
                 current_file_face2='',
                 gbc_file='',
                 name = 'lens_PO'):
        self.name = name
        self.freqList = freq
        self.lens = lens
        self.get_field = get_field
        self.method = method
        self.waist_radius = waist_radius
        self.po_points =po_points
        self.factor = factor
        self.spill_over = spill_over
        if coor_sys == '':
            self.coor_sys_name = coor_sys
            self.coor_sys =coor_sys
        else:
            self.coor_sys =coor_sys
            self.coor_sys_name = coor_sys.name
        if current_file_face1 == '':
            self.curr_file_face1 = self.name+'_face1.cur'
        else:
            self.curr_file_face1 =current_file_face1
        if current_file_face2 == '':
            self.curr_file_face2 = self.name+'_face2.cur'
        else:
            self.curr_file_face2 =current_file_face2
        self.curr_file_face2 = current_file_face2
        self.gbc_file = gbc_file

        self.Str = write_lens_po(self.name,self.freqList.name, 
                                 self.lens.name, 
                                 get_field=self.get_field,
                                 method=self.method, waist_radius=self.waist_radius,
                                 po_points=self.po_points, factor=self.factor,
                                 spill_over=self.spill_over,coord_sys=self.coor_sys_name, 
                                 current_file_face1=self.curr_file_face1, 
                                 current_file_face2=self.curr_file_face2,
                                 gbc_file=self.gbc_file)
        
class aperture_po():
    def __init__(self,
                 freq,
                 aperture,
                 method='po_plus_ptd',
                 po_points=[0,0], ptd_points= [[-1,0]],
                 factor=[0,0],
                 spill_over='off',
                 ray_output='none',
                 coor_sys='', 
                 file_name='',
                 name = 'Aperture_PO'):
        self.name = name
        self.freqList = freq
        self.scatter = aperture
        self.method = method
        self.po_points =po_points
        self.ptd_points = ptd_points
        self.factor = factor
        self.spill_over = spill_over
        self.ray_output = ray_output
        self.file_name = file_name
        if coor_sys == '':
            self.coor_sys_name = coor_sys
            self.coor_sys =coor_sys
        else:
            self.coor_sys =coor_sys
            self.coor_sys_name = coor_sys.name
        if file_name == '':
            self.file_name= self.name+'.cur'
        else:
            self.file_name= self.file_name
        self.Str = write_Aperture_PO(self.name,self.freqList.name, 
                                 self.scatter.name, 
                                 method=self.method,
                                 po_points=self.po_points, 
                                 ptd_points = self.ptd_points,
                                 factor=self.factor,
                                 spill_over=self.spill_over,
                                 ray_output = self.ray_output,
                                 coord_sys=self.coor_sys_name,
                                 current_filename =self.file_name)

class Spherical_grid():
    def __init__(self,
                 coor_sys,
                 u_range,v_range,
                 u0,v0,Nu,Nv,
                 grid_type= 'uv',
                 Truncation = 'rectangular',
                 e_h = 'e_field',
                 polarisation = 'linear',
                 near_far='near',
                 near_dist=100,
                 filename='',
                 name='spher_grid'):
        self.name = name
        self.coor_sys = coor_sys
        self.Beam_center = [u0,v0]
        self.Beam_size =[u_range,v_range]
        self.Beam_sampleN = [Nu,Nv]
        self.near_far = near_far
        self.near_dist = near_dist
        self.outputfile = filename
        self.grid_type = grid_type
        self.Truncation =Truncation
        self.e_h = e_h
        self.polarisation = polarisation


        self.Str = write_spherical_grid(self.name,
                                        self.coor_sys.name,
                                        self.Beam_size[0],self.Beam_size[1],
                                        self.Beam_center[0],self.Beam_center[1],
                                        self.Beam_sampleN[0],self.Beam_sampleN[1],
                                        grid_type = self.grid_type,
                                        Truncation = self.Truncation,
                                        e_h = self.e_h,
                                        polarisation =self.polarisation,
                                        near_far=self.near_far,
                                        near_dist=self.near_dist,
                                        filename=self.outputfile)
#%%
class Spherical_cut():
    def __init__(self,
                 coor_sys,
                 theta_range = [-10,10,501],
                 phi_range=[0,90,3],
                 cut_type = 'polar',
                 e_h = 'e_field',
                 polarisation = 'linear',
                 near_far='near',
                 near_dist=100,
                 filename='',
                 name='spher_cut'):
        self.name = name
        self.coor_sys = coor_sys
        self.theta_range = theta_range
        self.phi_range = phi_range
        self.cut_type =cut_type
        self.e_h = e_h
        self.polarisation = polarisation
        self.near_far = near_far
        self.near_dist = near_dist
        self.outputfile = filename


        self.Str = write_spherical_cut(self.name,
                                       self.coor_sys.name,
                                       theta_range = self.theta_range,
                                       phi_range = self.phi_range,
                                       cut_type =self.cut_type,
                                       e_h = self.e_h,
                                       polarisation = self.polarisation,
                                       near_far=self.near_far,
                                       near_dist=self.near_dist,
                                       filename=self.outputfile)        

