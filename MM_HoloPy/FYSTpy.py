import numpy as np
import torch as T

#geometry packages
from coordinate import coord_sys
from surface import PolySurf
from rim import Table_rect_rim
from elements import Reflector
#electral packages
from Feedpy import Gaussianfeed
from Kirchhoffpy import kirchhoff

''' Function used to read the data of the CCAT-prime Telescope
    & Holographic configurations.
'''
def read_FYSTholo(filename):
    pass
