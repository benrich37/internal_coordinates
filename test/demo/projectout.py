from ase.io import read,write
import numpy as np
import os
import sys
mainpath = os.path.join(os.path.dirname(__file__), '../../')
sys.path.append(mainpath)
from funcs import wilson, ic_gen, misc

proponal = read('../toy_data/proponal.xyz')
proponal2 = read('../toy_data/proponal2.xyz')

#disp_vec = proponal2.get_positions() - proponal.get_positions()
bonds_dict = misc.bonds_dict(proponal)
print(bonds_dict)