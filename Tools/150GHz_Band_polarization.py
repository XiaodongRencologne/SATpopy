from pro_grasp import SATpy_v2
import os
import numpy as np

"""output folder"""
folder = 'output/sim_bandpass/polarization_study/'
if not os.path.exists(folder):
    os.mkdir(folder)
beam_folder = '"D:/Xiaodong/SO/SAT_simple/SAT_simulations/sim_bandpass/150GHz_Band/cut/'

"""Input data"""
# Input 1 Rx position in focal plane:

Rx_x_list = [0,20,60,100,160,180]

polar_angle_list = [0]#np.arange(0,50,5)
freq_list = [90,150,220]
for Rx_x in Rx_x_list:
    Rx = [Rx_x,0,0]
    Folder = folder +'Rx_'+str(Rx[0])+'y_'+str(Rx[1])+'/'
    if not os.path.exists(Folder):
        os.mkdir(Folder)
    # Input 2: frequency list and beam profile files
    
    print(Folder)

    for item in freq_list:
        if not os.path.exists(Folder+str(item)+'/'):
            os.mkdir(Folder+str(item)+'/')
        print(Folder)
        print(Folder + str(item) +'/')
        f_eff = 569.56
        Ax = np.arctan(Rx[0]/f_eff)
        Ay = np.arctan(Rx[1]/f_eff)
        for polar_angle in polar_angle_list:
            if not os.path.exists(Folder+str(item)+'/'+str(polar_angle)+'degree/'):
                os.mkdir(Folder+str(item)+'/'+str(polar_angle)+'degree/')
            Model = SATpy_v2.SAT_v2(item,
                                    polar_angle,
                                    beam_folder + str(item)+'_GHz.cut"',
                                    Rx,
                                    method = {'lens':'go_plus_po',
                                            'screen': 'po_plus_ptd'},
                                    grids = {'type': 'uv',
                                            'center':[-Ax,Ay],
                                            'range': [0.4,0.4],
                                            'Points': [501,501]},
                                    cuts  = {'type': 'polar',
                                            'theta_center':[0.0,0.0],
                                            'theta_range': [-1,1,11],
                                            'phi_range':[0,90,3]},
                                    outputfolder = Folder+str(item)+'/'+str(polar_angle)+'degree/',
                                    srf_folder = 'D:/Xiaodong/SO/SAT_simple/SAT_simulations/sim_bandpass/150GHz_Band/srf/')
            Model._write_tor_tci()