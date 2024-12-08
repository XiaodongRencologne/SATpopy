from pathlib import Path

from run import run_command
from SATpy import SAT_v1
optics = SAT_v1('150GHz',[0,0,0])
folder_parent = 'output\\'
folders= ['test1\\','test2\\','test3\\']
folders = [folder_parent+item for item in folders ]

for item in folders:
    Path(item).mkdir(parents=True,exist_ok=True)
    optics =SAT_v1('150GHz',[0,0,0],outputfolder=item)
    optics._write_tor_tci()
    del(optics)
    #command ='copy '+str(Path.cwd())+'\\batch.gxp'+' '+str(Path.cwd())+ '\\' +item
    command = 'copy batch.gxp ' + item
    print(command)
    run_command(["cmd", "/c", "copy", "batch.gxp", item])
    run_command('bash test_sh.sh')