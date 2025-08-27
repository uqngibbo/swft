"""
Wrapper layer for testing the fitting of a scrf sim to CFD data.

@author: Nick Gibbons
"""

from numpy import sqrt, array, absolute, arange
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import subprocess
from test_friction import read_solution_file
import yaml
import matplotlib.pyplot as plt

f_ref = 1e-3
H_ref = 1e7

TEMPLATE = """
gas_file_name: "gm.lua"
reaction_file_name: "rr.lua"
xi: [0.0, 1.0]
Ai: [0.007853981633974483, 0.007853981633974483]
xf: [0.0, 1.0]
f:  [{}, {}]
dt: 5e-7
Hdot: {}
T0: 361.0
p0: 968.0
v0: 3623.0
Y0: {{"air": 1.0}}
calc_derivatives: True
"""

def objective(x):
    f0= x[0]*f_ref
    f1= x[1]*f_ref
    H = x[2]*H_ref
    with open('temp.yaml', 'w') as fp:
        fp.write(TEMPLATE.format(f0, f1, H))

    cmd = "scrf temp.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    current = read_solution_file('temp.bin')

    n = current[var].size
    interp = interp1d(target['x'], target[var])
    interped_v = interp(current['x'])
    current_v = current[var]

    L22= ((current_v-interped_v)**2).sum()/n
    print("Called obj with f=[{},{}] H={} and output={}".format(f0/f_ref, f1/f_ref, H/H_ref, L22))
    return L22

def jacobian(x):
    f0= x[0]*f_ref
    f1= x[1]*f_ref
    H = x[2]*H_ref
    with open('temp.yaml', 'w') as fp:
        fp.write(TEMPLATE.format(f0, f1, H))

    cmd = "scrf temp.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    current = read_solution_file('temp.bin')

    n = current[var].size
    interp = interp1d(target['x'], target[var])
    interped_v = interp(current['x'])
    current_v = current[var]

    fderivs0 = read_solution_file("fderivs_0000-temp.bin")
    fderivs1 = read_solution_file("fderivs_0001-temp.bin")
    Jf0 = 2.0/n*((current_v-interped_v)*fderivs0[var]*f_ref).sum()
    Jf1 = 2.0/n*((current_v-interped_v)*fderivs1[var]*f_ref).sum()

    Hderivs = read_solution_file("Hderivs-temp.bin")
    Hderivs_v = Hderivs[var]
    JH = 2.0/n*((current_v-interped_v)*Hderivs_v*H_ref).sum()

    print("       jac with f=[{},{}] H={} and output={},{},{}".format(f0/f_ref,f1/f_ref,H/H_ref,Jf0,Jf1,JH))
    return array([Jf0, Jf1, JH])

#var = 'v'
var = 'M'
target = read_solution_file('combined.bin')

#f0 = 1.0
#H = 5.0
##for f in [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010]: 
#for f1 in arange(1.0, 10.0, 1.0): 
#    out1 = objective(array([f0+1e-3, f1, H]))
#    out0 = objective(array([f0, f1, H]))
#    deriv = (out1-out0)/1e-3
#    jac = jacobian(array([f0, f1, H]))
#
#    print("Out1 {} out0 {}".format(out1, out0))
#    print("deriv {}".format(deriv))
#    print("jac   {}".format(jac[0]))

f0= 2.0
f1= 4.0
H = 7.0
x0 = array([f0, f1, H])
print("Starting with f=[{:5.5f},{:5.5f}] H={:5.5e}".format(f0, f1, H))
#res = minimize(objective, x0, bounds=(array([0.0, 100.0]),array([0.0, 100.0]),array([0.0, 20.0])), method='SLSQP', jac='2-point', options={'eps': 1e-3,'finite_diff_rel_step':1e-3}, tol=1e-8)
res = minimize(objective, x0, bounds=(array([0.0, 100.0]),array([0.0, 100.0]),array([0.0, 20.0])), method='SLSQP', jac=jacobian, tol=1e-8, options={'disp':True})

print("Done... found f0,f1,H={})".format(res.x))
print(res)

current = read_solution_file('temp.bin')
fig = plt.figure(figsize=(14,4))
axes0,axes1,axes2 = fig.subplots(1,3)

markersize = 3.0
axes0.plot(target['x'], target['v'], color='goldenrod', label="Optimiser Sol.", linewidth=2.0)
axes0.plot(current['x'][::10], current['v'][::10], color='black', linestyle='None', marker='o', markerfacecolor="None", markersize=markersize, label="Reference Sol.")
axes0.set_xlabel("x (m)")
axes0.set_title("Velocity (m/s)")
axes0.grid()
axes0.legend(framealpha=1.0)

axes1.plot(target['x'], target['T'], 'r-', linewidth=2.0)
axes1.plot(current['x'][::10], current['T'][::10],color='black', linestyle='None', marker='o',  markerfacecolor="None", markersize=markersize)
axes1.set_xlabel("x (m)")
axes1.set_title("Temperature (K)")
axes1.grid()

axes2.plot(target['x'], target['M'], color='cyan', linewidth=2.0)
axes2.plot(current['x'][::10], current['M'][::10],color='black', linestyle='None', marker='o', markerfacecolor="None", markersize=markersize)
axes2.set_xlabel("x (m)")
axes2.set_title("Mach Number (M)")
axes2.grid()

plt.tight_layout()
plt.show()
