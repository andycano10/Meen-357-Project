import subfunctions as sb
import matplotlib.pyplot as plt
from scipy.interpolate import *
from scipy.integrate import *
from statistics import mean

experiment,end_event = sb.experiment1()
rover, planet = sb.define_rover_1()


rover= sb.simulate_rover(rover, planet, experiment, end_event)
fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(6,8))
fig.suptitle('Rover Trajectory Graphs')

ax1.plot(rover['telemetry']['time'],rover['telemetry']['position'])
ax1.set_title('Position vs Time')
ax1.set_ylabel('Position (m)')
ax1.set_xlabel('Time (s)')

ax2.plot(rover['telemetry']['time'],rover['telemetry']['velocity'])
ax2.set_title('Velocity vs Time')
ax2.set_ylabel('Velocity (m/s)')
ax1.set_xlabel('Time (s)')


ax3.plot(rover['telemetry']['time'],rover['telemetry']['power'])
ax3.set_title('Power vs Time')
ax3.set_ylabel('Power (W)')
ax3.set_xlabel('Time (s)')

plt.tight_layout()
plt.show()
