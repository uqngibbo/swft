"""
Test code for swft. This is the heat addition only test.

@author: Nick
"""

from numpy import array, zeros, interp, frombuffer, array, concatenate, log, linspace, isclose
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import subprocess
from libswft import *

# Rayleigh flow with heat addition
def Tt_on_Ttstar(M, k):
    return (k+1)*M**2*(2+(k-1)*M**2)/(1+k*M**2)**2

def T_on_Tstar(M, k):
    return (k+1)**2*M**2/(1+k*M**2)**2

def p_on_pstar(M, k):
    return (k+1)**2*M**2/(1+k*M**2)**2

def rho_on_rhostar(M, k):
    return (1+k*M**2)/((k+1)*M**2)

def v_on_vstar(M, k):
    return ((k+1)*M**2)/(1+k*M**2)

def test_runswft():
    cmd = "swft ../examples/heat_addition.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

def test_output():
    data = read_solution_file('heat_addition.bin')
    assert isclose(data['M'][-1], 6.432, 1e-4)
    assert isclose(data['v'][-1], 3589.635, 1e-4)
    assert isclose(data['rho'][-1], 9.426e-3, 1e-4)
    assert isclose(data['p'][-1], 2097.011, 1e-4)
    assert isclose(data['T'][-1], 774.844, 1e-4)

def test_cleanup():
    cmd = "rm heat_addition.bin"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

if __name__=='__main__':
    data = read_solution_file('heat_addition.bin')
    gamma = data['gamma'][0]
    M0 = data['M'][0]

    pstar = data['p'][0]/p_on_pstar(M0, gamma)
    Tstar = data['T'][0]/T_on_Tstar(M0, gamma)
    rhostar = data['rho'][0]/rho_on_rhostar(M0, gamma)
    vstar = data['v'][0]/v_on_vstar(M0, gamma)

    print("pstar: ", pstar)
    print("Tstar: ", Tstar)
    print("rhostar: ", rhostar)

    every = 20
    Ms = data['M'][::every]
    ref = {}
    ref['x'] = data['x'][::every]
    ref['p'] = pstar*p_on_pstar(Ms, gamma)
    ref['T'] = Tstar*T_on_Tstar(Ms, gamma)
    ref['rho'] = rhostar*rho_on_rhostar(Ms, gamma)
    ref['v'] = vstar*v_on_vstar(Ms, gamma)
    ref['M'] = Ms


    fig = plt.figure(figsize=(10,8))
    axes = fig.subplots(2,2)

    axes[0,0].plot(data['x'], data['T'], 'r-')
    axes[0,0].plot(ref['x'], ref['T'], 'r.')
    axes[0,0].set_xlabel('x (m)')
    axes[0,0].set_ylabel('Temperature (K)')
    axes[0,0].grid()

    axes[0,1].plot(data['x'], data['rho'], 'g-')
    axes[0,1].plot(ref['x'], ref['rho'], 'g.')
    axes[0,1].set_xlabel('x (m)')
    axes[0,1].set_ylabel('density (kg/m3)')
    axes[0,1].grid()

    axes[1,0].plot(data['x'], data['v'], 'k-')
    axes[1,0].plot(ref['x'], ref['v'], 'k.')
    axes[1,0].set_xlabel('x (m)')
    axes[1,0].set_ylabel('Velocity (m/s)')
    axes[1,0].grid()

    axes[1,1].plot(data['x'], data['p'], 'b-')
    axes[1,1].plot(ref['x'], ref['p'], 'b.')
    axes[1,1].set_xlabel('x (m)')
    axes[1,1].set_ylabel('Pressure (Pa)')
    axes[1,1].grid()

    fig.suptitle("Compressible flow with heat addition")
    plt.tight_layout()
    plt.show()
