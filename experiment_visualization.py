import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import *
from scipy.integrate import *
import numpy as np


import numpy as np
def experiment1():
    experiment = {'time_range' : np.array([0,20000]),
    'initial_conditions' : np.array([0.3025,0]),
    'alpha_dist' : np.array([0, 100, 200, 300, 400, 500, 600, \
    700, 800, 900, 1000]),
    'alpha_deg' : np.array([11.509, 2.032, 7.182, 2.478, \
    5.511, 10.981, 5.601, -0.184, \
    0.714, 4.151, 4.042]),
    'Crr' : 0.1}
    # Below are default values for example only:
    end_event = {'max_distance' : 1000,
    'max_time' : 10000,
    'min_velocity' : 0.01}
    return experiment, end_event

experiment, end_event = experiment1()
alpha_dist = experiment['alpha_dist']
alpha_deg = experiment['alpha_deg']


alpha_fun = interp1d(alpha_dist, alpha_deg, kind = 'cubic', fill_value= 'extrapolate')

x_val = np.linspace(0,1000, 100)
angle = alpha_fun(x_val)

plt.plot(x_val, alpha_fun(x_val))
plt.plot(alpha_dist, alpha_deg, marker = '*', linestyle = 'None')
plt.title("Terrain angle vs. Position")
plt.xlabel('Position (m)')
plt.ylabel(r'Terrain Angle ($\theta$)')
plt.show()