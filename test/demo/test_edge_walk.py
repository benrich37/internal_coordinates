from ase.io import read,write
import numpy as np
import os
import sys
mainpath = os.path.join(os.path.dirname(__file__), '../../')
sys.path.append(mainpath)
from funcs import wilson, ic_gen, misc

proponal = read('../toy_data/proponal.xyz')
bonds_ref = misc.bonds_dict(proponal)
edges = misc.get_edges(bonds_ref)
print(edges)