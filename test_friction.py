"""
Test code for edict, currently conduction only.

@author: Nick
"""

from numpy import array, zeros, interp, frombuffer, array, concatenate, log, linspace
import struct
from io import BytesIO, BufferedWriter
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# Steady constant volume flow with wall friction

def f_Lstar_on_D(M, gamma):
    A = (1.0-M**2)/gamma/M**2
    B = (gamma+1)/2.0/gamma
    C = (gamma+1)*M**2/(2+(gamma-1)*M**2)
    print("Called with M={} f={}".format(M, A+B*log(C)))
    return A + B*log(C)

def p_on_pstar(M, gamma):
    return 1.0/M*((gamma+1)/(2+(gamma-1)*M**2))**0.5

def rho_on_rhostar(M, gamma):
    return 1.0/M*((2+(gamma-1)*M**2)/(gamma+1))**0.5

def T_on_Tstar(M, gamma):
    return (gamma+1)/(2+(gamma-1)*M**2)

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

f = 0.005
D = 0.05*2
L = 1.0
dLs = linspace(0.0, L, 20)
fLD = f_Lstar_on_D(M0, gamma)

Ms = []
for dL in dLs:
    constant = f*dL/D
    function = lambda Mout : fLD - constant - f_Lstar_on_D(Mout, gamma)
    Me = brentq(function, M0, 1.1)
    print("dL: {} Me: {}".format(dL, Me))
    Ms.append(Me)

Ms = array(Ms)
pstar = data['p'][0]/p_on_pstar(M0, gamma)
Tstar = data['T'][0]/T_on_Tstar(M0, gamma)
rhostar = data['rho'][0]/rho_on_rhostar(M0, gamma)

print("pstar: ", pstar)
print("Tstar: ", Tstar)
print("rhostar: ", rhostar)

ref = {}
ref['x'] = dLs
ref['p'] = pstar*p_on_pstar(Ms, gamma)
ref['T'] = Tstar*T_on_Tstar(Ms, gamma)
ref['rho'] = rhostar*rho_on_rhostar(Ms, gamma)
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

fig.suptitle("Compressible flow with friction (f={:3.3f})".format(f))
plt.tight_layout()
plt.show()
