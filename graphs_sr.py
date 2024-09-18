# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 20:30:27 2024

@author: catal
"""

import numpy as np
import matplotlib.pyplot as plt
from subfunctions import *

#plot 1: w_out (omega/Ng) speed  vs. t_out =t/Ng
r = get_gear_ratio(speed_reducer)


w_values = np.linspace(0,3.8,100)
w = []
for i in range(len(w_values)):
    w.append(w_values[i]) 
omega = np.array(w)

wout = []
for i in range(len(omega)):
    new_val = omega[i]/r
    wout.append(new_val)
    
tau = tau_dcmotor(omega, motor)
t = np.array(tau)
tau_out = []
for i in range(len(t)):
    new = t[i]*r
    tau_out.append(new)

#plt.plot(tau_out, wout)
 
#GRAPH 2:t_out vs. Pout [Nm] (use torque on the x-axis)
P = tau*omega
plt.plot(tau_dcmotor(omega, motor), P)
plt.show()



#GRAPH 3: motor power [W] vs. motor shaft speed [rad/s] (use speed on the x-axis)
plt.plot(wout, P)
plt.show()


