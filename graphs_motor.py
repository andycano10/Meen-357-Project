import numpy as np
import matplotlib.pyplot as plt

from subfunctions import *

#GRAPH 1: motor shaft speed [rad/s] vs. motor shaft torque [Nm] (use torque on the x-axis)

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
ax1.set_title('speed vs torque')
ax1.set_title('Speed vs Torque')
ax1.set_xlabel('Motor Shaft Torque [Nm]')
ax1.set_ylabel('Motor Speed [rad/s]')

ax2.plot(tau_dcmotor(omega, motor),P)
ax2.set_title('Power vs Torque')
ax2.set_xlabel('Motor Shaft Torque [Nm]')
ax2.set_ylabel('Motor Power [W]')


ax3.plot(w, P)
ax3.set_title('Motor Power vs Motor Speed')
ax3.set_xlabel('Motor Shaft Speed [rad/s]')
ax3.set_ylabel('Motor Power [W]')

plt.tight_layout()
plt.show()