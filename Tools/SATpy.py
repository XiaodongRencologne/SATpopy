#%%
import numpy as np

from GOElement import coor_sys, simple_lens, Aperture_screen, rim, global_coord
from EOElement import frequencyList, GaussBeam, GaussBeam_Near, Elliptical_Beam
from EOElement import lens_PO,aperture_po
from EOElement import Spherical_grid,Spherical_cut

from Command import get_current, get_field
SILICON = 3.36
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

#analysis_method_lens = 'go_plus_po'
#analysis_method_ap_po = 'po_plus_ptd'
lens_diameter1 = 44.8 # cm
lens_diameter2 = 44.8 # cm
lens_diameter3 = 44.8 # cm

class SAT_v1():
    def __init__(self, 
                 freq,
                 Rx_position,
                 method = {'lens':'go_plus_po',
                           'screen': 'po_plus_ptd'},
                 grids = {'type': 'uv',
                          'center':[0.0,0.0],
                          'range': [40,40],
                          'Points': [1001,1001]},
                 cuts  = {'type': 'polar',
                          'theta_center':[0.0,0.0],
                          'theta_range': [-90,90,4001],
                          'phi_range':[0,180,13]},
                 outputfolder=''):
        self.eff_focal_length = 569.56 #mm
        self.grids = grids
        self.cuts = cuts
        self.method =method
        self.outputfolder = outputfolder

        self.ap_list = []
        self.lens_list = []

        self.commit = {'Geometry': '',
                       'Feed': '',
                       'Frequency':'',
                       'Output':'',
                       'Methods':'',
                       'flow':''}

        self._create_coord(Rx_position=Rx_position)
        self._create_lens()
        self._create_rim()
        self._create_ap()
        # electrical objects
        self._create_input(freq)
        self._create_output()
        # methods
        self._create_Analysis()
        self._create_commands()
        #self._write_tor_tci()
       
    def _create_coord(self, Rx_position = [0,0,0], roatation = 0):
        ## 1. define coordinate systems
        self.coor_ref = coor_sys([0,0,0],[0,0,0],ref_coor = global_coord,name='coor_feed_ref',)
        self.coor_feed_offset = coor_sys(Rx_position,[np.pi,0,0],ref_coor = self.coor_ref,name='coor_feed_offset',)
        self.coor_feed_rot = coor_sys([0,0,0],[0,0,0],ref_coor = self.coor_feed_offset,name='coor_feed_rot',)
        self.coor_feed = coor_sys([0,0,0],[0,0,0],ref_coor = self.coor_feed_rot,name='coor_feed',)

        self.coor_lens3 = coor_sys([0,0,-L_lens3_ref*10],[0,0,0],ref_coor = self.coor_ref,name='coor_lens3',)
        self.coor_lens2 = coor_sys([0,0,-L_lens2_ref*10],[0,0,0],ref_coor = self.coor_ref,name='coor_lens2',)
        #coor_filter1 = coor_sys([0,0,-],[0,0,0],ref_coor = coor_lens2,name='coor_filter1')
        self.coor_lens1 = coor_sys([0,0,-L_lens1_ref*10],[0,0,0],ref_coor = self.coor_ref,name='coor_lens1',)
        self.coor_Lyot = coor_sys([0,0,-L_Ly_ref *10],[0,0,0],ref_coor = self.coor_ref,name='coor_Lyot',)
        #coor_filter2 = coor_sys([0,0,-804.9559555174151],[0,0,0],ref_coor = coor_ref,name='coor_filter',)
        #coor_IR = coor_sys([0,0,-816.9559555174151],[0,0,0],ref_coor = coor_ref,name='coor_IR',)
        self.coor_vw = coor_sys([0,0,-L_vw_ref*10],[0,0,0],ref_coor = self.coor_ref, name = 'coor_vw')
        self.coor_boresight_ref = coor_sys([0,0,0],[np.pi,0,0],ref_coor = self.coor_ref,name='coor_boresight_ref')
        Ax = -Rx_position[0]/self.eff_focal_length
        Ay = Rx_position[1]/self.eff_focal_length
        print(Ax/np.pi*180,Ay/np.pi*180)
        self.coor_cut = coor_sys([0,0,0],[Ay,Ax,0],ref_coor = self.coor_boresight_ref,name='coor_cut')
        self.coor_list = [global_coord, self.coor_ref, 
                            self.coor_feed_offset, self.coor_feed_rot, self.coor_feed, 
                            self.coor_lens1, self.coor_lens2, self.coor_lens3,
                            self.coor_Lyot,
                            self.coor_vw,
                            self.coor_boresight_ref, 
                            self.coor_cut]
        
    def _create_lens(self):
        ### 2. define lenses
        self.lens1 = simple_lens(self.coor_lens1, lens_diameter1, 
                                SILICON, loss_tangent = 0, 
                                r1 = '0 cm', r2 = '0 cm', bs1 = 0, bs2 = 0,
                                thickness = '4.34991 cm',#4.349908221542306 cm',
                                surf_f1 = '../srf/lens1_f1.rsf',
                                surf_f2 = '../srf/lens1_f2.rsf',
                                lengthUnit = 'cm',
                                name='lens1')

        self.lens2 = simple_lens(self.coor_lens2, lens_diameter2, 
                                SILICON, loss_tangent = 0, 
                                r1 = '0 cm', r2 = '0 cm', bs1 = 0, bs2 = 0,
                                thickness = '4.69671 cm', #'4.696706712699847 cm',
                                surf_f1 = '../srf/lens2_f1.rsf',
                                surf_f2 = '../srf/lens2_f2.rsf',
                                lengthUnit = 'cm',
                                name='lens2')

        self.lens3 = simple_lens(self.coor_lens3, lens_diameter3, 
                                SILICON, loss_tangent = 0, 
                                r1 = '0 cm', r2 = '0 cm', bs1 = 0, bs2 = 0,
                                thickness = '2.96556 cm' ,#'2.965564711384346 cm',
                                surf_f1 = '../srf/lens3_f1.rsf',
                                surf_f2 = '../srf/lens3_f2.rsf',
                                lengthUnit = 'cm',
                                name='lens3')        
        self.lens_list = [self.lens1,
                          self.lens2,
                          self.lens3]
        
    def _create_rim(self):
        self.rim_Lyot = rim([0,0], [210,210], Type = 'elliptical_rim',name='rim_Lyot')
        self.rim_vw = rim([0,0], [2.876618694117068E+002,2.876618694117068E+002], Type = 'elliptical_rim',name='rim_vw')
        self.rim_list = [self.rim_Lyot,
                         self.rim_vw]

    def _create_ap(self):
        self.Lyot = Aperture_screen(self.coor_Lyot, self.rim_Lyot,infinity_shadow = 'on', name='Lyot')
        self.VW = Aperture_screen(self.coor_vw, self.rim_vw,infinity_shadow = 'on', name='vw')
        self.ap_list=[self.Lyot,
                      self.VW]
        
    def _create_input(self,freq):
        self.freq_list = frequencyList([Feed_list[freq]['freq']], name ='freq_list')
        self.Feed_ellip = Elliptical_Beam(self.freq_list, self.coor_feed,
                                          Feed_list[freq]['Ellip_Taper'],
                                          [Feed_list[freq]['T_angle_x'],Feed_list[freq]['T_angle_y']],
                                          polarisation='linear',
                                          polarisation_angle=90,
                                          far_forced = 'off',
                                          factor = [0,0],
                                          frequency_index_for_plot = 1,
                                          name = 'Gaussian_Elliptical_Beam')
        
        self.Feed_Gaussian = GaussBeam(self.freq_list,self.coor_feed,
                                       Feed_list[freq]['T_angle_x'],
                                       Feed_list[freq]['Ellip_Taper'],
                                       polarisation='linear_y',
                                       name='Gauss_circle')
        
        self.input_list = [self.freq_list,
                           self.Feed_ellip,
                           self.Feed_Gaussian]
        
    def _create_Analysis(self):
        ### define PO analysis object
        self.PO_lens1 = lens_PO(self.freq_list,
                                self.lens1,get_field='lens_in_screen',
                                method=self.method['lens'], waist_radius=0,
                                po_points={'face1':[0,0],'face2': [0,0]}, factor=[0,0],
                                spill_over='on',coor_sys='', 
                                current_file_face1='', 
                                current_file_face2='',
                                gbc_file='',
                                name = 'lens1_PO')

        self.PO_lens2 =lens_PO(self.freq_list,
                        self.lens2,get_field='lens_in_screen',
                        method=self.method['lens'], waist_radius=0,
                        po_points={'face1':[0,0],'face2': [0,0]}, factor=[0,0],
                        spill_over='on',coor_sys='', 
                        current_file_face1='', 
                        current_file_face2='',
                        gbc_file='',
                        name = 'lens2_PO')

        self.PO_lens3 =lens_PO(self.freq_list,
                        self.lens3,get_field='lens_in_screen',
                        method=self.method['lens'], waist_radius=0,
                        po_points={'face1':[0,0],'face2': [0,0]}, factor=[0,0],
                        spill_over='on',coor_sys='', 
                        current_file_face1='', 
                        current_file_face2='',
                        gbc_file='',
                        name = 'lens3_PO')

        self.POA_Lyot = aperture_po(self.freq_list,
                            self.Lyot,
                            method=self.method['screen'],
                            po_points=[0,0], ptd_points= [[-1,0]],
                            factor=[0,0],
                            spill_over='on',
                            ray_output='none',
                            coor_sys='', 
                            file_name='',
                            name = 'POA_Lyot'
                            )
        
        self.POA_VW = aperture_po(self.freq_list,
                       self.VW,
                       method=self.method['screen'],
                       po_points=[0,0], ptd_points= [[-1,0]],
                       factor=[0,0],
                       spill_over='on',
                       ray_output='none',
                       coor_sys='', 
                       file_name='',
                       name = 'POA_VW'
                       )
        
        self.method_list = [self.PO_lens1,self.PO_lens2,self.PO_lens3,self.POA_Lyot,self.POA_VW] 
        pass
    def _create_output(self):
        self.Beam_grid = Spherical_grid(self.coor_cut,
                                        self.grids['range'][0],self.grids['range'][0],
                                        self.grids['center'][0], self.grids['center'][1],
                                        self.grids['Points'][0], self.grids['Points'][1],
                                        grid_type = self.grids['type'],
                                        Truncation = 'rectangular',
                                        e_h = 'e_field',
                                        polarisation = 'linear',
                                        near_far='far',
                                        filename='',
                                        name='Beam_grid')
        self.Beam_cut = Spherical_cut(self.coor_cut,
                                      theta_range = self.cuts['theta_range'],
                                      phi_range = self.cuts['phi_range'],
                                      cut_type = self.cuts['type'],
                                      e_h = 'e_field',
                                      polarisation = 'linear',
                                      near_far='far',
                                      filename='',
                                      name='Beam_cut')
        self.output_list = [self.Beam_cut,self.Beam_grid]
    
    def _create_commands(self):
        # create commands flow
        get_lens3_cur = get_current(self.PO_lens3,[self.Feed_ellip],
                                    accuracy= -80,
                                    auto_convergence=True,
                                    convergence_on_scatterer = [self.lens2])

        get_lens2_cur = get_current(self.PO_lens2,[self.PO_lens3],
                                    accuracy= -80,
                                    auto_convergence=True,
                                    convergence_on_scatterer = [self.lens1])

        get_lens1_cur = get_current(self.PO_lens1,[self.PO_lens2],
                                    accuracy= -80,
                                    auto_convergence=True,
                                    convergence_on_scatterer = [self.Lyot])

        get_Lyot_cur = get_current(self.POA_Lyot,[self.PO_lens1],
                                    accuracy= -80,
                                    auto_convergence=True,
                                    convergence_on_scatterer= [self.VW])

        get_VW_cur = get_current(self.POA_VW,[self.POA_Lyot],
                                accuracy= -80,
                                auto_convergence=True,
                                convergence_on_output_grid = [self.Beam_grid])


        get_field_grid = get_field(self.Beam_grid, [self.POA_VW])
        get_field_cut = get_field(self.Beam_cut, [self.POA_VW])

        self.command_list = [get_lens3_cur,get_lens2_cur,get_lens1_cur,
                             get_Lyot_cur, get_VW_cur,
                             get_field_grid,get_field_cut]
        
    def _write_tor_tci(self,):
        with open(self.outputfolder+'sat_optics.tor','w') as f:
            # write frequency
            #f.writelines(self.freq_list.Str)
            self.commit['Frequency']+=self.freq_list.name

            # write coordinates
            for item in self.coor_list:
                f.writelines(item.Str)

            # write lens
            for item in self.lens_list:
                f.writelines(item.Str)
                self.commit['Geometry']+=item.name+','
                
            # write rim
            for item in self.rim_list:
                f.writelines(item.Str)

            # write aperture screen
            for item in self.ap_list:
                f.writelines(item.Str)
                self.commit['Geometry']+=item.name+','
            
            # input Gaussian beams
            for item in self.input_list:
                f.writelines(item.Str)
                self.commit['Feed']+=item.name +','
            for item in self.output_list:
                f.writelines(item.Str)
                self.commit['Output']+=item.name +','
            
            # Analysis methods
            for item in self.method_list:
                f.writelines(item.Str)
                self.commit['Methods']+=item.name+','

        with open(self.outputfolder+'sat_optics.tci','w') as f:
            for item in self.command_list:
                f.writelines(item.Str)
                #self.commit['flow'] +=item.name
            pass
# %%
