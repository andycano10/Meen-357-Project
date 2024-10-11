import subfunctions as sb
import matplotlib.pyplot as plt


experiment,end_event = sb.experiment1()
end_event['max_distance'] = 1000
end_event['max_time'] = 10000
end_event['min_velocity'] = 0.01

rover = sb.simulate_rover(sb.rover, sb.planet, experiment, end_event)
fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(6,8))
fig.suptitle('Rover Trajectory Graphs')

ax1.plot(rover['telemetry']['position'],rover['telemetry']['time'])
ax1.set_title('Position vs Time')
ax1.set_xlabel('Position [m]')
ax1.set_ylabel('Time [s]')

ax2.plot(rover['telemetry']['velocity'],rover['telemetry']['time'])
ax2.set_title('Velocity vs Time')
ax2.set_xlabel('Velocity [m/s]')
ax1.set_ylabel('Time [s]')


ax3.plot(rover['telemetry']['power'],rover['telemetry']['time'])
ax3.set_title('Power vs Time')
ax3.set_xlabel('Power [W]')
ax1.set_ylabel('Time [s]')

plt.tight_layout()
plt.show()
