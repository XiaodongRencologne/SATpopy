#%%
import subprocess
import argparse
import time
# %%
#print(s[1])

def run_command(command):
    #import shlex
    #process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
    process = subprocess.Popen(command,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE,
                               text=True)
    for line in process.stdout:
        print(line)
    while True:
        output = process.stdout.readline()
        if process.poll() is not None:
            print(f"Process finished with exit code: {process.poll()}")
            break
        time.sleep(1)
    print('command completed!!')


def run_grasp(command):
    #import shlex
    #process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
    ex_cmd = '/cygdrive/c/Program\ Files/TICRA/TICRA-Tools-20.1.2/bin/ticra-tools.exe'
    process = subprocess.Popen([ex_cmd,command],
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE,
                               text=True)
    for line in process.stdout:
        print(line)
    while True:
        output = process.stdout.readline()
        if process.poll() is not None:
            print(f"Process finished with exit code: {process.poll()}")
            break
        time.sleep(1)
    print('command completed!!')
#%%

#run_command('python test.py')
#%%
#outputfolder=''
#ticra-tools.exe
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", type=str, default='')
    args = parser.parse_args()

    #opt1_value = args.opt1
    #opt2_value = args.opt2
    run_command('python test.py')
"""
