from pathlib import Path

from run import run_command
from One_Lens import Lens_v1

freq = '90GHz'
folder_parent = 'output\\Sim_Lens1\\'+freq+'\\'
method = {'lens':'po_plus_po',
          'screen': 'po_plus_ptd'}
methods = ['go_plus_po',
           'po_plus_po']

for item in methods:
    filefolder = folder_parent+item+'\\'
    Path(filefolder).mkdir(parents=True,exist_ok=True)
    method['lens'] = item
    optics =Lens_v1(freq,[0,0,0],method = method ,outputfolder=filefolder)
    optics._write_tor_tci()
    del(optics)
    #command ='copy '+str(Path.cwd())+'\\batch.gxp'+' '+str(Path.cwd())+ '\\' +item
    #command = 'copy batch.gxp ' + filefolder
    #print(command)
    #run_command(["cmd", "/c", "copy", "batch.gxp", filefolder])
    #run_command('bash test_sh.sh')