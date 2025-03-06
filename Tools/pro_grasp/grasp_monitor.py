import sys
try:
  import psutil
except ImportError:
   print("Cannot import psutil module - this is needed for this application.")
   print("Exiting...")
   sys.exit();
style="margin:0px auto"

import time
import re
import argparse as ap
import datetime

if __name__ == '__main__':

    parser = ap.ArgumentParser(formatter_class=\
                                         ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', action='store', dest='fname', type=str,
                        help='Filename for output',
                        default='grasp_monitor.log')


    parser.add_argument('-d', action='store', dest='fdir', type=str,
                        help='Directory for output log file',
                        default='../log/')

    parser.add_argument('-s', action='store', dest='sleeptime', type=int,
                        help='Number seconds to sleep between logging',
                        default=5)

    args = parser.parse_args()

    cnt = 0
    fid = open(args.fdir+args.fname, 'w+')

    while True:

        print(cnt)
        ts = time.time()

        #string = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        string = '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now() +
            datetime.timedelta(hours=+1))

        string += '\n'

        #perc = psutil.cpu_times_percent(interval=0.5, percpu=False)
        mem = psutil.virtual_memory();

        percs = psutil.cpu_percent(interval=0.5, percpu=False)
        string += "CPU User:     {:5.2f} %\n".format(percs)
        string += "Memory Used:     {:d} M\n".format(int(mem.used / 1024 / 1024))
        string += "Memory Percent:  {:5.2f}%\n".format(mem.percent)
        print(string)

        fid.write(string)
        fid.flush()
        time.sleep(args.sleeptime)

        cnt+=1

    fid.close()



