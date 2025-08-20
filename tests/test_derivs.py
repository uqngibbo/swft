"""
Test code for edict, currently conduction only.

@author: Nick
"""

from numpy import array, zeros, interp, frombuffer, array, concatenate, log, linspace, isclose
import struct
from io import BytesIO, BufferedWriter
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import subprocess
from test_friction import read_solution_file
import yaml

def read_cfg(filename):
    with open(filename) as fp:
        f = yaml.safe_load(fp)
    cfg = {}
    for k,v in f.items(): cfg[k] = v
    return cfg

name0 = 'friction'
cfg0  = read_cfg("{}.yaml".format(name0))
data0= read_solution_file("{}.bin".format(name0))

name1 = 'perturb_friction'
cfg1  = read_cfg("{}.yaml".format(name1))
data1= read_solution_file("{}.bin".format(name1))

keys = ['p', 'v', 'rho']
dUdf_fd = {key:(data1[key]-data0[key])/(cfg1['f']-cfg0['f']) for  key in keys}

dUdf = read_solution_file("derivs-{}.bin".format(name0))

for key in keys:
    print(key)
    print(("dUdf_fd"+" {:12.12e}"*len(dUdf_fd[key])).format(*dUdf_fd[key]))
    print(("dUdf   "+" {:12.12e}"*len(dUdf[key])).format(*dUdf[key]))
