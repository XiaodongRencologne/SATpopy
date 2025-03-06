import numpy as np

'''
1. read cut files
'''
def read_cut(filename):
    cuts = {}
    with open(filename,'r') as f:
        while True:
            line = f.readline()  # Read a line
            if not line:  # Stop if end of file is reached
                break
            elif line == 'Tabulated feed data\n':
                    inf = f.readline().split()
                    phi = str(int(float(inf[3])))
                    cuts[phi] = {}
                    theta = np.linspace(float(inf[0]), float(inf[0]) + float(inf[1]) * (float(inf[2])-1), int(inf[2]))
                    cuts[phi]['theta'] = theta
                    data = np.genfromtxt(f,max_rows = int(inf[2]))
                    cuts[phi]['E_co'] = data[:,0]
                    cuts[phi]['E_cx'] = data[:,1]
                    '''
                    for _ in range(int(inf[2])):
                         if not f.readline():
                              break
                    '''
    return cuts