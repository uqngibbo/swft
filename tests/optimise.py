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
from pylab import plot,show

f_ref = 1e-3
H_ref = 1e7

TEMPLATE = """
gas_file_name: "gm.lua"
reaction_file_name: "rr.lua"
L: 1.0
rs: 0.05
rf: 0.05
dt: 5e-7
f: {}
Hdot: {}
T0: 361.0
p0: 968.0
v0: 3623.0
Y0: {{"air": 1.0}}
calc_derivatives: True
"""

def objective(x):
    f = x[0]*f_ref
    H = x[1]*H_ref
    with open('temp.yaml', 'w') as fp:
        fp.write(TEMPLATE.format(f, H))

    cmd = "scrf temp.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    current = read_solution_file('temp.bin')
    #ref = absolute(target['v']).max()
    #ref = 1.0
    # Not sure about this reference value

    n = current['v'].size
    interp = interp1d(target['x'], target['v'])
    interped_v = interp(current['x'])
    current_v = current['v']

    L22= ((current_v-interped_v)**2).sum()/n
    print("Called obj with f={} H={} and output={}".format(f/f_ref, H/H_ref, L22))
    return L22

def jacobian(x):
    f = x[0]*f_ref
    H = x[1]*H_ref
    with open('temp.yaml', 'w') as fp:
        fp.write(TEMPLATE.format(f, H))

    cmd = "scrf temp.yaml"
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)
    assert proc.returncode == 0, "Failed cmd: "+cmd

    current = read_solution_file('temp.bin')
    #ref = absolute(target['v']).max()
    #ref = 1.0

    n = current['v'].size
    interp = interp1d(target['x'], target['v'])
    interped_v = interp(current['x'])
    current_v = current['v']

    fderivs = read_solution_file("fderivs-temp.bin")
    fderivs_v = fderivs['v']
    Jf = 2.0/n*((current_v-interped_v)*fderivs_v*f_ref).sum()

    Hderivs = read_solution_file("Hderivs-temp.bin")
    Hderivs_v = Hderivs['v']
    JH = 2.0/n*((current_v-interped_v)*Hderivs_v*H_ref).sum()

    print("       jac with f={} H={} and output={},{}".format(f/f_ref,H/H_ref,Jf,JH))
    return array([Jf, JH])

target = read_solution_file('combined.bin')

#f = 0.0
##for f in [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010]: 
#for H in arange(1.0, 10.0, 1.0): 
#    out1 = objective(array([f, H+1e-3]))
#    out0 = objective(array([f, H]))
#    deriv = (out1-out0)/1e-3
#    jac = jacobian(array([f, H]))
#
#    print("Out1 {} out0 {}".format(out1, out0))
#    print("deriv {}".format(deriv))
#    print("jac   {}".format(jac))

f = 1.0
H = 7.0
x0 = array([f, H])
print("Starting with f={:5.5f} H={:5.5e}".format(f, H))
#res = minimize(objective, x0, bounds=(array([1e-4, 0.1]),), method='SLSQP', jac=jacobian, options={'eps':   10e-6,'finite_diff_rel_step':1e-1}, tol=1e-6)
res = minimize(objective, x0, bounds=(array([0.0, 100.0]),array([0.0, 20.0])), method='SLSQP', jac=jacobian, tol=1e-9, options={'disp':True})

print("Done... found f,H={})".format(res.x))
print(res)

current = read_solution_file('temp.bin')
plot(current['x'][::10], current['v'][::10], 'k.')
plot(target['x'], target['v'], 'r-')
show()
