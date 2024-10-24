# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 20:30:27 2024

@author: catal
"""

import numpy as np
import matplotlib.pyplot as plt
from subfunctions import *

r = get_gear_ratio(speed_reducer)


w_values = np.linspace(0,3.8,100)
w = []
for i in range(len(w_values)):
    w.append(w_values[i]) 
omega = np.array(w)
wout = omega/r

    
tau = tau_dcmotor(omega, motor)
t = np.array(tau)
tau_out = t*r



P = tau_out*wout

fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(6,8))
fig.suptitle('Graphs for the Speed Reducer')

#plot 1: w_out (omega/Ng) speed  vs. t_out =t/Ng
ax1.plot(tau_out, wout)
ax1.set_title('Speed vs Torque')
ax1.set_xlabel('Reducer Torque [Nm]')
ax1.set_ylabel('Reducer Speed [rad/s]')

#GRAPH 2:t_out vs. Pout [Nm] (use torque on the x-axis)
ax2.plot(tau_out,P)
ax2.set_title('Power vs Torque')
ax2.set_xlabel('Reducer Torque [Nm]')
ax2.set_ylabel('Reducer Power [W]')


#GRAPH 3: motor power [W] vs. motor shaft speed [rad/s] (use speed on the x-axis)
ax3.plot(wout, P)
ax3.set_title('Power vs Speed')
ax3.set_xlabel('Reducer Speed [rad/s]')
ax3.set_ylabel('Reducer Power [W]')

<<<<<<< Updated upstream
plt.tight_layout()
=======
plt.tight_layout()
plt.show()

plt.tight_layout()


Pmax_index = np.argmax(P)
speed = wout[50]
print(max(P), Pmax_index, speed)



>>>>>>> Stashed changes
