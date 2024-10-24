#%%
import numpy as np
import W2GRASP_GO, W2GRASP_EO
from W2GRASP_EO import write_frequency_list, write_Gauss_beam, write_Gauss_beam_near,write_Gauss_Ellip_Beam
from W2GRASP_EO import write_lens_po
from W2GRASP_EO import write_spherical_grid
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
        self.Str = write_Gauss_beam_near(self.name,self.freqList,
                                         self.coor_sys,
                                         self.b_radius, self.phase_radius,
                                         polarisation = self.polarisation,
                                         factor = self.factor,
                                         frequency_index_for_plot = self.freq_index)

class Elliptical_Beam():
    def __init__(self,
                 freqList,
                 coor_sys,
                 Taper, T_angle, 
                 polarisation, polarisation_angle,
                 far_forced,
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
                                          self.freqList,
                                          self.coor_sys,
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
        else:
            self.coor_sys =coor_sys
            self.coor_sys_name = coor_sys.name
        self.coor_sys = coor_sys
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
        

class Spherical_grid():
    def __init__(self,
                 coor_sys,
                 u_range,v_range,
                 u0,v0,Nu,Nv,
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


        self.Str = write_spherical_grid(self.name,
                                        self.coor_sys.name,
                                        self.Beam_size[0],self.Beam_size[1],
                                        self.Beam_center[0],self.Beam_center[1],
                                        self.Beam_sampleN[0],self.Beam_sampleN[1],
                                        near_far=self.near_far,
                                        near_dist=self.near_dist,
                                        filename=self.outputfile)
#%%


