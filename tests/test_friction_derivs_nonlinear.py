"""
Test code for edict, currently conduction only.

@author: Nick
"""

from numpy import array, zeros, interp, frombuffer, array, concatenate, log, linspace, isclose, sqrt
import struct
from io import BytesIO, BufferedWriter
import matplotlib.pyplot as plt
import subprocess
from test_friction import read_solution_file
import yaml
from scipy.interpolate import interp1d

TEMPLATE = """
gas_file_name: "gm.lua"
reaction_file_name: "rr.lua"
xi: [0.0, 1.0]
Ai: [0.007853981633974483, 0.007853981633974483]
xf: {}
f:  {}
dt: 5e-7
Hdot: 0.0
T0: 361.0
p0: 968.0
v0: 3623.0
Y0: {{"air": 1.0}}
calc_derivatives: {}
"""

def read_cfg(filename):
    with open(filename) as fp:
        f = yaml.safe_load(fp)
    cfg = {}
    for k,v in f.items(): cfg[k] = v
    return cfg

def read_derivs(name0, name1):
    cfg0  = read_cfg("{}.yaml".format(name0))
    data0= read_solution_file("{}.bin".format(name0))
    fi0= interp1d(cfg0['xf'], cfg0['f'])
    data0['f'] = fi0(data0['x'])

    cfg1  = read_cfg("{}.yaml".format(name1))
    data1= read_solution_file("{}.bin".format(name1))
    fi1= interp1d(cfg1['xf'], cfg1['f'])
    data1['f'] = fi1(data1['x'])

    #keys = ['p', 'v', 'rho', 'M', 'T']
    #dUdf_fd = {key:(data1[key]-data0[key])/(f1-f0) for  key in keys}

    dUdf = read_solution_file("fderivs-{}.bin".format(name0))
    return data0, data1, dUdf

def get_L2_norms(dUdf_fd, dUdf):
    n = dUdf_fd['rho'].size
    L2 = {}
    L2['rho'] = sqrt(((dUdf_fd['rho']-dUdf['rho'])**2).sum()/n)
    L2['p']   = sqrt(((dUdf_fd['p']-dUdf['p'])**2).sum()/n)
    L2['v']   = sqrt(((dUdf_fd['v']-dUdf['v'])**2).sum()/n)
    return L2

def test_runscrf():
    xf= [0.0, 1.0]
    f = [0.005, 0.002]
    f2= [0.005001, 0.002001]
    with open('temp.yaml', 'w') as fp:
        fp.write(TEMPLATE.format(xf, f, "true"))

    cmd = "scrf temp.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    with open('temp_perturbed.yaml', 'w') as fp:
        fp.write(TEMPLATE.format(xf, f2, "false"))

    cmd = "scrf temp_perturbed.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

def test_output():
    data0, dUdf_fd, dUdf = read_derivs()
    L2 = get_L2_norms(dUdf_fd, dUdf)

    assert isclose(L2['rho'], 4.6954e-6, 1e-4)
    assert isclose(L2['p'], 11.184, 1e-4)
    assert isclose(L2['v'], 0.35706, 1e-4)

def test_cleanup():
    cmd = "rm temp.bin temp_perturbed.bin Hderivs-temp.bin fderivs-temp.bin"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd


if __name__=='__main__':
    #xf= [0.0, 1.0]
    #f = [0.005, 0.002]
    #f2= [0.005001, 0.002001]
    #with open('temp.yaml', 'w') as fp:
    #    fp.write(TEMPLATE.format(xf, f, "true"))

    #cmd = "scrf temp.yaml"
    #proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    #assert proc.returncode == 0, "Failed cmd: "+cmd

    #with open('temp_perturbed.yaml', 'w') as fp:
    #    fp.write(TEMPLATE.format(xf, f2, "false"))

    #cmd = "scrf temp_perturbed.yaml"
    #proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    #assert proc.returncode == 0, "Failed cmd: "+cmd

    data0, data1, dUdf = read_derivs('baseline', 'perturbed')

    df = 1e-6
    keys = ['p', 'v', 'rho', 'M', 'T']
    dUdf_fd = {key:(data1[key]-data0[key])/(df) for key in keys}
    print("Density derivative with respect to f0: ")
    print(dUdf_fd['rho'])
    print("velocity derivative with respect to f0: ")
    print(dUdf_fd['v'])
    print("pressure derivative with respect to f0: ")
    print(dUdf_fd['p'])

    #L2 = get_L2_norms(dUdf_fd, dUdf)
    #print("rho L2 norm: ", L2['rho'])
    #print("p   L2 norm: ", L2['p'])
    #print("v   L2 norm: ", L2['v'])

    fig = plt.figure(figsize=(16,4))
    axes = fig.subplots(1,5)

    axes[0].plot(data0['x'], dUdf_fd['rho'], linestyle='--', marker='.', linewidth=2.0, color='red')
    axes[0].plot(data0['x'], dUdf['rho'], linestyle='-', marker='.', linewidth=1.0, color='maroon')
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

    axes[3].plot(data0['x'], dUdf_fd['M'], linestyle='--', linewidth=2.0, color='grey')
    axes[3].plot(data0['x'], dUdf['M'], linestyle='-', linewidth=1.0, color='black')
    axes[3].set_xlabel('x (m)')
    axes[3].set_ylabel('dMdf')
    axes[3].grid()

    axes[4].plot(data0['x'], dUdf_fd['T'], linestyle='--', linewidth=2.0, color='orange')
    axes[4].plot(data0['x'], dUdf['T'], linestyle='-', linewidth=1.0, color='red')
    axes[4].set_xlabel('x (m)')
    axes[4].set_ylabel('dTdf')
    axes[4].grid()

    fig.suptitle("Method of Total Derivatives Check")
    plt.tight_layout()
    plt.show()

