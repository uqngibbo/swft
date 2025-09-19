"""
Test code for swft, checking basic Q1D area change

@author: Nick
"""

from numpy import array, zeros, interp, frombuffer, array, concatenate, isclose
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import subprocess
from libswft import *

# Validate using steady-isentropic flow formulas from Frank White, Ch 9

def A_on_Astar(M, gamma):
    a = 1.0 + 0.5*(gamma-1)*M**2
    b = 0.5*(gamma+1)
    c = 0.5*(gamma+1)/(gamma-1)
    return 1.0/M*(a/b)**c

def A_on_A0(M0, M, gamma):
    top = 1/M* ((1+0.5*(gamma-1)*M**2) /(0.5*(gamma+1)))**(0.5*(gamma+1)/(gamma-1))
    bot = 1/M0*((1+0.5*(gamma-1)*M0**2)/(0.5*(gamma+1)))**(0.5*(gamma+1)/(gamma-1))
    return top/bot

def T_on_T0(M0, M, gamma):
    top = (1.0 + (gamma-1)/2.0*M**2)**-1
    bot = (1.0 + (gamma-1)/2.0*M0**2)**-1
    return top/bot

def p_on_p0(M0, M, gamma):
    top = (1.0 + (gamma-1)/2.0*M**2)**(-gamma/(gamma-1))
    bot = (1.0 + (gamma-1)/2.0*M0**2)**(-gamma/(gamma-1))
    return top/bot

def rho_on_rho0(M0, M, gamma):
    top = (1.0 + (gamma-1)/2.0*M**2)**(-1.0/(gamma-1))
    bot = (1.0 + (gamma-1)/2.0*M0**2)**(-1.0/(gamma-1))
    return top/bot

def test_runswft():
    cmd = "swft ../examples/area_change.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

def test_output():
    data = read_solution_file('area_change.bin')

    assert isclose(data['M'][-1], 12.7455, 1e-4)
    assert isclose(data['v'][-1], 3665.695, 1e-4)
    assert isclose(data['rho'][-1], 2.29823e-3, 1e-4)
    assert isclose(data['p'][-1], 135.786, 1e-4)
    assert isclose(data['T'][-1], 205.789, 1e-4)

def test_cleanup():
    cmd = "rm area_change.bin"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd
    
if __name__=='__main__':
    data = read_solution_file('area_change.bin')
    gamma = data['gamma'][0]
    M0 = data['M'][0]
    A0 = data['A'][0]
    Astar = A0/A_on_Astar(M0, gamma)

    ref = {}
    Mref = []
    every = 10
    for Ai in data['A'][::every]:
        function = lambda Mguess : Ai/Astar - A_on_Astar(Mguess, gamma)
        Me = brentq(function, 1.0, 16.0)
        Mref.append(Me)

    Mref = array(Mref)
    ref['M'] = Mref
    ref['x'] = data['x'][::every]
    ref['A'] = data['A'][0]*A_on_A0(M0, Mref, gamma)
    ref['T'] = data['T'][0]*T_on_T0(M0, Mref, gamma)
    ref['p'] = data['p'][0]*p_on_p0(M0, Mref, gamma)
    ref['rho'] = data['rho'][0]*rho_on_rho0(M0, Mref, gamma)
    
    print("L2 error M:   ", (((ref['M']-data['M'][::every])**2).sum()/Mref.size)**0.5)
    print("L2 error T:   ", (((ref['T']-data['T'][::every])**2).sum()/Mref.size)**0.5)
    print("L2 error p:   ", (((ref['p']-data['p'][::every])**2).sum()/Mref.size)**0.5)
    print("L2 error rho: ", (((ref['rho']-data['rho'][::every])**2).sum()/Mref.size)**0.5)

    print("data['M'][-1]", data['M'][-1])
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

    axes[1,0].plot(data['x'], data['M'], 'k-')
    axes[1,0].plot(ref['x'], ref['M'], 'k.')
    axes[1,0].set_xlabel('x (m)')
    axes[1,0].set_ylabel('Mach Number')
    axes[1,0].grid()

    axes[1,1].plot(data['x'], data['p'], 'b-')
    axes[1,1].plot(ref['x'], ref['p'], 'b.')
    axes[1,1].set_xlabel('x (m)')
    axes[1,1].set_ylabel('Pressure (Pa)')
    axes[1,1].grid()

    plt.tight_layout()
    plt.show()
