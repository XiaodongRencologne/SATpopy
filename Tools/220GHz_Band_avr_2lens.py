from pro_grasp import SATpy_v2_two_lens
import os
import numpy as np

"""output folder"""
folder = 'output/sim_bandpass/220GHz_Band/'
if not os.path.exists(folder):
    os.mkdir(folder)
beam_folder = '"D:/Xiaodong/SO/SAT_simple/SAT_simulations/sim_bandpass/220GHz_Band/cut/'

"""Input data"""
# Input 1 Rx position in focal plane:

Rx_x_list = np.arange(0,201,10)

for Rx_x in Rx_x_list:
    Rx = [Rx_x,0,0]
    Folder = folder +'Rx_'+str(Rx[0])+'y_'+str(Rx[1])+'/'
    if not os.path.exists(Folder):
        os.mkdir(Folder)
    # Input 2: frequency list and beam profile files
    freq_list = [220]#np.linspace(220,221,1,dtype = np.int32)
    print(Folder)

    for item in freq_list:
        if not os.path.exists(Folder+str(item)+'/'):
            os.mkdir(Folder+str(item)+'/')
        print(Folder)
        print(Folder + str(item) +'/')
        f_eff = 569.56
        Ax = np.arctan(Rx[0]/f_eff)*180/np.pi
        Ay = np.arctan(Rx[1]/f_eff)*180/np.pi
        Model = SATpy_v2_two_lens.SAT_v2(item,
                                    beam_folder + str(item)+'_GHz.cut"',
                                    Rx,
                                    grids = {'type': 'EloverAz',
                                            'center':[Ax,Ay],
                                            'range': [10,10],
                                            'Points': [2001,2001]},
                                    cuts  = {'type': 'polar',
                                            'theta_center':[0.0,0.0],
                                            'theta_range': [-5,5,4001],
                                            'phi_range':[0,90,3]},
                                    outputfolder = Folder+str(item)+'/',
                                    srf_folder = 'D:/Xiaodong/SO/SAT_simple/SAT_simulations/sim_bandpass/220GHz_Band/srf/')
        Model._write_tor_tci()