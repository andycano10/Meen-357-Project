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


# plt.figure()
# plt.subplot(3,1,1)
# plt.plot(tau_dcmotor(omega, motor), w_values)
# plt.xlabel("Motor Shaft Torque (Nm)")
# plt.ylabel("motor shaft \n speed [rad/s]")
# plt.title("Sped vs Torque")



#GRAPH 2:motor power [W] vs. motor shaft torque [Nm] (use torque on the x-axis)
tau = tau_dcmotor(omega, motor)
# t = np.array(tau)
P = tau * omega

# plt.subplot(3,1,2)
# plt.plot(tau_dcmotor(omega, motor), P)
# plt.xlabel("Motor Shaft Torque (Nm)")
# plt.ylabel("Motor Power [W]")
# plt.ylim(0,200)

# # #GRAPH 3: motor power [W] vs. motor shaft speed [rad/s] (use speed on the x-axis)
# plt.subplot(3,1,3)
# plt.plot(w, P)
# plt.xlabel("motor shaft speed [rad/s]")
# plt.ylabel("motor power [W]")
# plt.ylim(0,200)
# plt.tight_layout()
# plt.show()
fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(6,8))
fig.suptitle('test')
ax1.plot(tau_dcmotor(omega, motor), omega)
ax2.plot(tau_dcmotor(omega, motor),P)
ax3.plot(w, P)
plt.tight_layout()
plt.show()
