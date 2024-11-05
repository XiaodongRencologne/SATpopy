#%%
import numpy as np

from GOElement import coor_sys, simple_lens, Aperture_screen, rim, global_coord
from EOElement import frequencyList, GaussBeam, GaussBeam_Near, Elliptical_Beam
from EOElement import lens_PO,aperture_po
from EOElement import Spherical_grid,Spherical_cut

from Command import get_current, get_field

L_lensFp_3   = 7.177590111674096
L_lens3_2    = 15.586806616226909
L_lens2_1    = 57.632802785493645
L_lens1_Lyot = 1.162050628144469
L_Ly_vw      = 22.7114

L_lens1_ref = L_lensFp_3 + L_lens3_2 + L_lens2_1
L_lens2_ref = L_lensFp_3 + L_lens3_2
L_lens3_ref = L_lensFp_3 
L_Ly_ref = L_lens1_ref + L_lens1_Lyot
L_vw_ref = L_Ly_ref + L_Ly_vw

Freq_list = ['90GHz', '150GHz' , '220GHz', '280GHz']
Feed_list = {'90GHz': {'freq': 90,
                       'Ellip_Taper': -2.1714724,
                       'T_angle_x' : 15.4270683,
                       'T_angle_y' : 20.2798893,
                       'Gauss_Taper':  -20*np.log10(np.exp(1)),
                       'Gauss_Tangle': 30.8,
                       'beam_radius': '1.972 mm',
                       'phase_front_radius': '0 mm'},
            '150GHz': {'freq':150,
                       'Ellip_Taper': -2.1714724,
                       'T_angle_x' : 10.1161095,
                       'T_angle_y' : 11.8095882,
                       'Gauss_Taper': -20*np.log10(np.exp(1)),
                       'Gauss_Tangle':20.232219,
                       'beam_radius': '1.802 mm',
                       'phase_front_radius': '0 mm'},
            '220GHz': {'freq':220,
                       'Ellip_Taper': -2.1714724,
                       'T_angle_x' : 7.060976,
                       'T_angle_y' : 7.9273612,
                       'Gauss_Taper': -20*np.log10(np.exp(1)),
                       'Gauss_Tangle': 14.121952,
                       'beam_radius': '1.760 mm',
                       'phase_front_radius': '0 mm'},
            '280GHz': {'freq':280,
                       'Ellip_Taper': -2.1714724,
                       'T_angle_x' : 5.7901323,
                       'T_angle_y' : 6.3342926,
                       'Gauss_Taper': -20*np.log10(np.exp(1)),
                       'Gauss_Tangle': 11.58,
                       'beam_radius': '1.686 mm',
                       'phase_front_radius': '0 mm'}
            }

grid_type = 'El&Az'
Center_angle = [0.0, 0.0]
x_range = 40
y_range = 40
Nx = 4001
Ny = 4001

cut_type = 'polar'
theta_range = [-90,90,8001]
phi_range = [0,180,13]

analysis_method_lens = 'go_plus_po'
analysis_method_ap_po = 'po_plus_ptd'
lens_diameter1 = 44.8 # cm
lens_diameter2 = 44.8 # cm
lens_diameter3 = 44.8 # cm

class SAT_v1():
    def __init__(self,
                 freq_list, 
                 feed_list,
                 Rx_position_list,
                 method = {'lens':'go_plus_po',
                           'screen': 'po_plus_ptd'},
                 grids = {'type': 'El&Az',
                          'center':[[0.0,0.0]],
                          'range': 40,
                          'Points': 4001},
                 cuts  = {'type': 'polar',
                          'theta_center':[[0.0,0.0]],
                          'theta_range': 90,
                          'theta_Points': 8001,
                          'phi_range':[0,180],
                          'cut_number':13},
                 outputfolder=''):
        
        self.commit = {'Geometry': 'Only 3 lenses',
                       'Feed': 'elliptical Gaussian Beam, force to far field',
                       }
        self._create_coord()
        self._create_lens()
        self._create_ap()
        self._create_rim()
        self._create_EO()


    def _create_coord(self):
        ## 1. define coordinate systems
        self.coor_ref = coor_sys([0,0,0],[0,0,0],ref_coor = global_coord,name='coor_feed_ref',)
        self.coor_feed_offset = coor_sys([0,0,0],[np.pi,0,0],ref_coor = self.coor_ref,name='coor_feed_offset',)
        self.coor_feed_rot = coor_sys([0,0,0],[0,0,0],ref_coor = self.coor_feed_offset,name='coor_feed_rot',)
        self.coor_feed = coor_sys([0,0,0],[0,0,0],ref_coor = self.coor_feed_rot,name='coor_feed',)

        self.coor_lens3 = coor_sys([0,0,-L_lensFp_3*10],[0,0,0],ref_coor = self.coor_ref,name='coor_lens3',)
        self.coor_lens2 = coor_sys([0,0,-L_lens3_2*10-L_lensFp_3*10],[0,0,0],ref_coor = self.coor_ref,name='coor_lens2',)
        #coor_filter1 = coor_sys([0,0,-],[0,0,0],ref_coor = coor_lens2,name='coor_filter1')
        self.coor_lens1 = coor_sys([0,0,-L_lens2_1*10-L_lens3_2*10-L_lensFp_3*10],[0,0,0],ref_coor = self.coor_ref,name='coor_lens1',)
        self.coor_Lyot = coor_sys([0,0,-L_lens1_Lyot*10-L_lens2_1*10-L_lens3_2*10-L_lensFp_3*10],[0,0,0],ref_coor = self.coor_ref,name='coor_Lyot',)
        #coor_filter2 = coor_sys([0,0,-804.9559555174151],[0,0,0],ref_coor = coor_ref,name='coor_filter',)
        #coor_IR = coor_sys([0,0,-816.9559555174151],[0,0,0],ref_coor = coor_ref,name='coor_IR',)
        self.coor_vw = coor_sys([0,0,-L_Ly_vw*10-L_lens1_Lyot*10-L_lens2_1*10-L_lens3_2*10-L_lensFp_3*10],[0,0,0],ref_coor = self.coor_ref, name = 'coor_vw')
        self.coor_cut = coor_sys([0,0,0],[np.pi,0,0],ref_coor = self.coor_ref,name='coor_cut',)
        self.coor_list=[global_coord, self.coor_ref, 
                        self.coor_feed_offset, self.coor_feed_rot, self.coor_feed, 
                        self.coor_lens1, self.coor_lens2, 
                        self.coor_lens3,
                        self.coor_Lyot,
                        self.coor_vw, 
                        self.coor_cut]
        
    def _create_lens(self):
        pass
        
# %%
test = SAT_v1(Freq_list,Feed_list,[0])
# %%
test.coor_list
# %%
