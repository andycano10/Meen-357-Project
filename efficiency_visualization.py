#By submitting this assignment, I agree to the following:
#   "Aggies do not lie, cheat, or steal, or tolerate those who do."
#   "I have not given or received any unauthorized aid on this assignment."
#
# Name:        Jackson Walton, Andy Cano-Avila, Catalina Fossati
# Section:      505
# Assignment:   THE ASSIGNMENT NUMBER (e.g. Lab 1b-2)
# Date:         DAY MONTH 2024

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import *
import subfunctions as sb

rover, planet = sb.define_rover_1()
effcy_tau = rover['wheel_assembly']['motor']['effcy_tau']
effcy = rover['wheel_assembly']['motor']['effcy']

effcy_fun = interp1d(effcy_tau, effcy, kind = 'cubic')

x0 = np.linspace(0, 165, 100)
effcy_val = effcy_fun(x0)

plt.plot(x0, effcy_val)
plt.plot(effcy_tau, effcy, '*')

plt.xlabel("Torque (N*m)")
plt.ylabel("Efficiency (unitless)")
plt.title("Efficiency vs Torque")
plt.show()