# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 12:49:19 2024

@author: catal
"""
import numpy as np
import matplotlib.pyplot as plt

from subfunctions import *

#GRAPH 1: motor shaft speed [rad/s] vs. motor shaft torque [Nm] (use torque on the x-axis)
#linespace for y

w_values = np.linspace(0,3.75,100)
w = []
for i in range(len(w_values)):
    w.append(w_values[i]) 
omega = np.array(w)


#GRAPH 2:motor power [W] vs. motor shaft torque [Nm] (use torque on the x-axis)
tau = tau_dcmotor(omega, motor)

P = tau * omega


# # #GRAPH 3: motor power [W] vs. motor shaft speed [rad/s] (use speed on the x-axis)
fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(6,8))
fig.suptitle('Graphs for the Motor')
ax1.plot(tau_dcmotor(omega, motor), omega)
ax1.ylabel('Speed Reducer')
ax1.xlabel('Motor Torque')
ax1.Title('Speed vs. Torque')
ax2.plot(tau_dcmotor(omega, motor),P)
ax3.plot(w, P)
plt.tight_layout()
plt.show()
