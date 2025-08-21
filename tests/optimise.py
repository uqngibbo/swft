"""
Wrapper layer for testing the fitting of a scrf sim to CFD data.

@author: Nick Gibbons
"""

from numpy import sqrt, array, absolute
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import subprocess
from test_friction import read_solution_file
import yaml
from pylab import plot,show


TEMPLATE = """
gas_file_name: "gm.lua"
reaction_file_name: "rr.lua"
L: 1.0
rs: 0.05
rf: 0.05
dt: 5e-7
f: {}
Hdot: 0.0
T0: 361.0
p0: 968.0
v0: 3623.0
Y0: {{"air": 1.0}}
calc_derivatives: True
"""

def objective(x):
    f = x[0]
    with open('temp.yaml', 'w') as fp:
        fp.write(TEMPLATE.format(f))

    cmd = "scrf temp.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    current = read_solution_file('temp.bin')
    ref = absolute(target['v']).max()
    # Not sure about this reference value

    n = current['v'].size
    interp = interp1d(target['x'], target['v'])
    interped_v = interp(current['x'])/ref
    current_v = current['v']/ref

    L22= ((current_v-interped_v)**2).sum()/n
    #print("Called obj with f={} and output={}".format(f, L22))
    return L22

def jacobian(x):
    f = x[0]

    cmd = "scrf temp.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    current = read_solution_file('temp.bin')
    ref = absolute(target['v']).max()

    n = current['v'].size
    interp = interp1d(target['x'], target['v'])
    interped_v = interp(current['x'])/ref
    current_v = current['v']/ref

    derivs = read_solution_file("fderivs-temp.bin")
    derivs_v = derivs['v']/ref

    Jf = 2.0/n*((current_v-interped_v)*derivs_v).sum()
    #print("Called jac with f={} and output={}".format(f, Jf))
    return array([Jf])

target = read_solution_file('friction.bin')

#for f in [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010]: 
#    out1 = objective(array([f+1e-5]))
#    out0 = objective(array([f]))
#    deriv = (out1-out0)/1e-5
#    jac = jacobian(array([f]))
#
#    print("Out1 {} out0 {}".format(out1, out0))
#    print("deriv {}".format(deriv))
#    print("jac   {}".format(jac))

f = 0.020
x0 = array([f])
print("Starting with f={:5.5f}".format(f))
#res = minimize(objective, x0, bounds=(array([1e-4, 0.1]),), method='SLSQP', jac=jacobian, options={'eps':   10e-6,'finite_diff_rel_step':1e-1}, tol=1e-6)
res = minimize(objective, x0, bounds=(array([1e-4, 0.1]),), method='SLSQP', jac=jacobian, tol=1e-9)

print("Done... found f={:5.5f} (should be {:5.5f})".format(res.x[0], 0.005))
print(res)

current = read_solution_file('temp.bin')
plot(current['x'][::10], current['v'][::10], 'k.')
plot(target['x'], target['v'], 'r-')
show()
