"""
Collection of swft related code in python.

@author: Nick Gibbons
"""

from numpy import array, zeros, interp, frombuffer, array, concatenate, log, linspace, isclose
import struct
from io import BytesIO, BufferedWriter
import subprocess
import yaml
from sys import argv
import matplotlib.pyplot as plt

def read_config_file(filename):
    """
    It's nice to have a little yaml file in each simulation directory, that stores
    some information about what's happening with that sim. This routine tries to
    read that file, and if not, returns a default configuation dictionary to use
    for plotting things.
    """

    config = {}
    with open(filename) as fp:
        f = yaml.safe_load(fp)
    for k,v in f.items(): config[k] = v
    return config

def read_solution_file(filename):
    if not filename.endswith('.bin'):
        raise Exception("Wrong file extension for file {}".format(filename))

    with open(filename, 'rb') as fp:
        bytes = fp.read()

    stream = BytesIO(bytes)

    buff = stream.read(8*2)
    neq, N = struct.unpack("Q"*2, buff)

    buff = stream.read(8*N); x     = frombuffer(buff)
    buff = stream.read(8*N); p     = frombuffer(buff)
    buff = stream.read(8*N); T     = frombuffer(buff)
    buff = stream.read(8*N); rho   = frombuffer(buff)
    buff = stream.read(8*N); A     = frombuffer(buff)
    buff = stream.read(8*N); v     = frombuffer(buff)
    buff = stream.read(8*N); M     = frombuffer(buff)
    buff = stream.read(8*N); gamma = frombuffer(buff)

    data = {}
    data['neq'] = neq; data['N'] = N;
    data['x']     = x.copy()
    data['p']     = p.copy()
    data['T']     = T.copy()
    data['rho']   = rho.copy()
    data['A']     = A.copy()
    data['v']     = v.copy()
    data['M']     = M.copy()
    data['gamma'] = gamma.copy()

    return data

filename = argv[1]
data = read_solution_file(filename)

fig = plt.figure(figsize=(10,8))
axes = fig.subplots(2,2)

axes[0,0].plot(data['x'], data['T'], 'r-')
axes[0,0].set_xlabel('x (m)')
axes[0,0].set_ylabel('Temperature (K)')
axes[0,0].grid()

axes[0,1].plot(data['x'], data['rho'], 'g-')
axes[0,1].set_xlabel('x (m)')
axes[0,1].set_ylabel('density (kg/m3)')
axes[0,1].grid()

axes[1,0].plot(data['x'], data['M'], 'k-')
axes[1,0].set_xlabel('x (m)')
axes[1,0].set_ylabel('Mach Number')
axes[1,0].grid()

axes[1,1].plot(data['x'], data['p'], 'b-')
axes[1,1].set_xlabel('x (m)')
axes[1,1].set_ylabel('Pressure (Pa)')
axes[1,1].grid()

plt.tight_layout()
plt.show()
