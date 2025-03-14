from pro_grasp import SATpy_v2
import os
import numpy as np

"""output folder"""
folder = 'output/sim_bandpass/150GHz_Band/'
if not os.path.exists(folder):
    os.mkdir(folder)
beam_folder = '"../cut/'

"""Input data"""
# Input 1 Rx position in focal plane:
Rx = [0,0,0]
folder = folder +'Rx_'+str(Rx[0])+'y_'+str(Rx[1])+'/'
if not os.path.exists(folder):
    os.mkdir(folder)
# Input 2: frequency list and beam profile files
freq_list = np.linspace(132,168,37,dtype = np.int32)

for item in freq_list:
    if not os.path.exists(folder+str(item)+'/'):
        os.mkdir(folder+str(item)+'/')
    Model = SATpy_v2.SAT_v2(item,
                            beam_folder + str(item)+'_GHz.cut"',
                            Rx,
                            grids = {'type': 'EloverAz',
                                      'center':[0.0,0.0],
                                      'range': [10,10],
                                      'Points': [2001,2001]},
                            cuts  = {'type': 'polar',
                                      'theta_center':[0.0,0.0],
                                      'theta_range': [-5,5,4001],
                                      'phi_range':[0,90,3]},
                            outputfolder = folder+str(item)+'/')
    Model._write_tor_tci()