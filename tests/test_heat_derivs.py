"""
Test code for swft, checking heat transfer derivs.

@author: Nick
"""

from numpy import array, zeros, interp, frombuffer, array, concatenate, log, linspace, isclose, sqrt
import matplotlib.pyplot as plt
import subprocess
import yaml
from libswft import *

def read_cfg(filename):
    with open(filename) as fp:
        f = yaml.safe_load(fp)
    cfg = {}
    for k,v in f.items(): cfg[k] = v
    return cfg

def read_derivs():
    name0 = 'baseline'
    cfg0  = read_cfg("{}.yaml".format(name0))
    data0= read_solution_file("{}.bin".format(name0))

    name1 = 'perturbed'
    cfg1  = read_cfg("{}.yaml".format(name1))
    data1= read_solution_file("{}.bin".format(name1))

    keys = ['p', 'v', 'rho']
    dUdH_fd = {key:(data1[key]-data0[key])/(cfg1['Hdot']-cfg0['Hdot']) for  key in keys}

    dUdH = read_solution_file("Hderivs-{}.bin".format(name0))
    return data0, dUdH_fd, dUdH

def get_L2_norms(dUdH_fd, dUdH):
    n = dUdH_fd['rho'].size
    L2 = {}
    L2['rho'] = sqrt(((dUdH_fd['rho']-dUdH['rho'])**2).sum()/n)
    L2['p']   = sqrt(((dUdH_fd['p']  -dUdH['p'])**2).sum()/n)
    L2['v']   = sqrt(((dUdH_fd['v']  -dUdH['v'])**2).sum()/n)
    return L2

def test_runswft():
    baseline = read_config_file('../examples/heat_addition.yaml')
    baseline['calc_derivatives'] = "true"
    write_config_to_file(baseline, 'baseline.yaml')

    cmd = "swft baseline.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    perturbed = {key:val for key,val in baseline.items()}
    perturbed['calc_derivatives'] = "false"
    perturbed['Hdot'] *= 1.001
    write_config_to_file(perturbed, 'perturbed.yaml')

    cmd = "swft perturbed.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

def test_output():
    data0, dUdH_fd, dUdH = read_derivs()
    L2 = get_L2_norms(dUdH_fd, dUdH)

    assert isclose(L2['rho'], 2.984e-16, 1e-4)
    assert isclose(L2['p'], 2.132e-09, 1e-4)
    assert isclose(L2['v'], 6.290e-11, 1e-4)

def test_cleanup():
    cmd = "rm baseline.bin perturbed.bin"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    cmd = "rm Hderivs-baseline.bin fderivs_0000-baseline.bin fderivs_0001-baseline.bin CHderivs-baseline.bin"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    cmd = "rm baseline.yaml perturbed.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd


if __name__=='__main__':

    data0, dUdH_fd, dUdH = read_derivs()

    L2 = get_L2_norms(dUdH_fd, dUdH)
    print("rho L2 norm: ", L2['rho'])
    print("p   L2 norm: ", L2['p'])
    print("v   L2 norm: ", L2['v'])

    fig = plt.figure(figsize=(14,4))
    axes = fig.subplots(1,3)

    axes[0].plot(data0['x'], dUdH_fd['rho'], linestyle='--', linewidth=2.0, color='red')
    axes[0].plot(data0['x'], dUdH['rho'], linestyle='-', linewidth=1.0, color='maroon')
    axes[0].set_xlabel('x (m)')
    axes[0].set_ylabel('drhodf')
    axes[0].grid()

    axes[1].plot(data0['x'], dUdH_fd['p'], linestyle='--', linewidth=2.0, color='cyan')
    axes[1].plot(data0['x'], dUdH['p'], linestyle='-', linewidth=1.0, color='blue')
    axes[1].set_xlabel('x (m)')
    axes[1].set_ylabel('dpdf')
    axes[1].grid()

    axes[2].plot(data0['x'], dUdH_fd['v'], linestyle='--', linewidth=2.0, color='green')
    axes[2].plot(data0['x'], dUdH['v'], linestyle='-', linewidth=1.0, color='forestgreen')
    axes[2].set_xlabel('x (m)')
    axes[2].set_ylabel('dvdf')
    axes[2].grid()

    fig.suptitle("Method of Total Derivatives Check")
    plt.tight_layout()
    plt.show()

