"""
Wrapper layer for testing the fitting of a scrf sim to CFD data.

@author: Nick Gibbons
"""

from numpy import sqrt, array
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

    print("Called obj with f={}".format(f))
    cmd = "scrf temp.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    current = read_solution_file('temp.bin')
    derivs = read_solution_file("derivs-temp.bin")

    n = current['v'].size
    interp = interp1d(target['x'], target['v'])
    interped_v = interp(current['x'])

    L2 = sqrt(((current['v']-interped_v)**2).sum()/n)
    print("L2[v]: ", L2)
    return L2

f = 0.010
x0 = array([f])
target = read_solution_file('friction.bin')

res = minimize(objective, x0, bounds=(array([1e-4, 0.1]),), method='SLSQP', jac='2-point', options={'eps':   10e-6,'finite_diff_rel_step':1e-1}, tol=1e-6)

print("Done... result is:")
print(res)

current = read_solution_file('temp.bin')
plot(current['x'][::10], current['v'][::10], 'k.')
plot(target['x'], target['v'], 'r-')
show()
