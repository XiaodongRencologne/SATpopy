from pathlib import Path

from run import run_command
from SATpy import SAT_v1
#optics = SAT_v1('150GHz',[0,0,0])
freq = '150GHz'
folder_parent = 'output\\'+freq+'\\'+'PO_PO'+'\\'
method = {'lens':'po_plus_po',
          'screen': 'po_plus_ptd'}
Rxs = [[0,0,0],
       [50,0,0],
       [90,0,0],
       [120,0,0],
       [140,0,0],
       [150,0,0],
       [160,0,0],
       [170,0,0],
       [180,0,0]]

for item in Rxs:
    filefolder = folder_parent+'_x'+str(item[0])+'_y'+str(item[1])+'\\'
    Path(filefolder).mkdir(parents=True,exist_ok=True)
    optics =SAT_v1(freq,
                   item,
                   method = method,
                   outputfolder=filefolder)
    optics._write_tor_tci()
    del(optics)
    #command ='copy '+str(Path.cwd())+'\\batch.gxp'+' '+str(Path.cwd())+ '\\' +item
    #command = 'copy batch.gxp ' + filefolder
    #print(command)
    #run_command(["cmd", "/c", "copy", "batch.gxp", filefolder])
    #run_command('bash test_sh.sh')