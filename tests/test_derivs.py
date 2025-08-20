"""
Test code for edict, currently conduction only.

@author: Nick
"""

from numpy import array, zeros, interp, frombuffer, array, concatenate, log, linspace, isclose, sqrt
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

def read_derivs():
    name0 = 'friction'
    cfg0  = read_cfg("{}.yaml".format(name0))
    data0= read_solution_file("{}.bin".format(name0))

    name1 = 'perturb_friction'
    cfg1  = read_cfg("{}.yaml".format(name1))
    data1= read_solution_file("{}.bin".format(name1))

    keys = ['p', 'v', 'rho']
    dUdf_fd = {key:(data1[key]-data0[key])/(cfg1['f']-cfg0['f']) for  key in keys}

    dUdf = read_solution_file("derivs-{}.bin".format(name0))
    return data0, dUdf_fd, dUdf

def get_L2_norms(dUdf_fd, dUdf):
    n = dUdf_fd['rho'].size
    L2 = {}
    L2['rho'] = sqrt(((dUdf_fd['rho']-dUdf['rho'])**2).sum()/n)
    L2['p']   = sqrt(((dUdf_fd['p']-dUdf['p'])**2).sum()/n)
    L2['v']   = sqrt(((dUdf_fd['v']-dUdf['v'])**2).sum()/n)
    return L2

def test_runscrf():
    cmd = "scrf friction.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    cmd = "scrf perturb_friction.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

def test_output():
    data0, dUdf_fd, dUdf = read_derivs()
    L2 = get_L2_norms(dUdf_fd, dUdf)

    assert isclose(L2['rho'], 6.1493e-06, 1e-4)
    assert isclose(L2['p'], 14.398, 1e-4)
    assert isclose(L2['v'], 0.45423, 1e-4)

def test_cleanup():
    cmd = "rm friction.bin perturb_friction.bin derivs-friction.bin"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd


if __name__=='__main__':
    data0, dUdf_fd, dUdf = read_derivs()

    L2 = get_L2_norms(dUdf_fd, dUdf)
    print("rho L2 norm: ", L2['rho'])
    print("p   L2 norm: ", L2['p'])
    print("v   L2 norm: ", L2['v'])

    fig = plt.figure(figsize=(14,4))
    axes = fig.subplots(1,3)

    axes[0].plot(data0['x'], dUdf_fd['rho'], linestyle='--', linewidth=2.0, color='red')
    axes[0].plot(data0['x'], dUdf['rho'], linestyle='-', linewidth=1.0, color='maroon')
    axes[0].set_xlabel('x (m)')
    axes[0].set_ylabel('drhodf')
    axes[0].grid()

    axes[1].plot(data0['x'], dUdf_fd['p'], linestyle='--', linewidth=2.0, color='cyan')
    axes[1].plot(data0['x'], dUdf['p'], linestyle='-', linewidth=1.0, color='blue')
    axes[1].set_xlabel('x (m)')
    axes[1].set_ylabel('dpdf')
    axes[1].grid()

    axes[2].plot(data0['x'], dUdf_fd['v'], linestyle='--', linewidth=2.0, color='green')
    axes[2].plot(data0['x'], dUdf['v'], linestyle='-', linewidth=1.0, color='forestgreen')
    axes[2].set_xlabel('x (m)')
    axes[2].set_ylabel('dvdf')
    axes[2].grid()

    fig.suptitle("Method of Total Derivatives Check")
    plt.tight_layout()
    plt.show()

