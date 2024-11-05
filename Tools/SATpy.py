import numpy as np

from GOElement import coor_sys, simple_lens, Aperture_screen, rim, global_coord
from EOElement import frequencyList, GaussBeam, GaussBeam_Near, Elliptical_Beam
from EOElement import lens_PO,aperture_po
from EOElement import Spherical_grid,Spherical_cut

from Command import get_current, get_field

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

grid_type = 'uv'
Center_angle = [0.0, 0.0]
x_range = 20
y_range = 20

cut_type = 'polar'
theta_range = [-10,10,1001]
phi_range = [0,90,3]

class SAT_v1():
    def __init__(self, freq_list, 
                 Feed_list,
                 method = 'po',
                 outputfolder=''):
        
        self.commit = {'Geometry': 'Only 3 lenses',
                       'Feed': 'elliptical Gaussian Beam, force to far field',
                       }
        