"""
Test code for edict, currently conduction only.

@author: Nick
"""

from numpy import array, zeros, interp, frombuffer, array, concatenate
import struct
from io import BytesIO, BufferedWriter
import matplotlib.pyplot as plt

# Validate using steady-isentropic flow formulas from Frank White, Ch 9

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

def read_solution_file(filename):
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
    data['x']     = x
    data['p']     = p
    data['T']     = T
    data['rho']   = rho
    data['A']     = A
    data['v']     = v
    data['M']     = M
    data['gamma'] = gamma

    return data

data = read_solution_file('solution.bin')
gamma = data['gamma'][0]
M0 = data['M'][0]

ref = {}
every = 10
ref['M'] = data['M'][::every]
ref['x'] = data['x'][::every]
Mref = ref['M']
ref['A'] = data['A'][0]*A_on_A0(M0, Mref, gamma)
ref['T'] = data['T'][0]*T_on_T0(M0, Mref, gamma)
ref['p'] = data['p'][0]*p_on_p0(M0, Mref, gamma)
ref['rho'] = data['rho'][0]*rho_on_rho0(M0, Mref, gamma)


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

axes[1,0].plot(data['x'], data['A'], 'k-')
axes[1,0].plot(ref['x'], ref['A'], 'k.')
axes[1,0].set_xlabel('x (m)')
axes[1,0].set_ylabel('Area (m2)')
axes[1,0].grid()

axes[1,1].plot(data['x'], data['p'], 'b-')
axes[1,1].plot(ref['x'], ref['p'], 'b.')
axes[1,1].set_xlabel('x (m)')
axes[1,1].set_ylabel('Pressure (Pa)')
axes[1,1].grid()

plt.tight_layout()
plt.show()
