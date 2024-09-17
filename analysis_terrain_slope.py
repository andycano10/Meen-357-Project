import numpy as np
import subfunctions as sb
import matplotlib.pyplot as plt

crr = 0.2 # coefficient of resistance
slope_array_deg = np.linspace(-10, 35, 25) # gives 25 degree angles between -10 and 35
v_max = [] # list of max velocities at different angles
iter_num_list = [] # list of number of iterations
err_est_list = [] # list of error estimate for root
exitFlag_list = [] # list of exit flags

def bisection(x): #find roots for F_net (F=0, v=max)

    lb = 0
    ub = 3.8 # m/s
    err_max = 1e-3
    iter_max = 1000

    # y-values for lower and upper bound
    fun_l = sb.F_net(np.array([lb]), np.array([x]), sb.rover, sb.planet, crr)
    fun_u = sb.F_net(np.array([ub]), np.array([x]), sb.rover, sb.planet, crr)

    # initate
    done = False
    iter_num = 0
    err_est = None
    root = None

    # incase no root
    if fun_l * fun_u > 0:
        exitFlag = -1 # no root
        velocity = root
        return velocity, err_est, iter_num, exitFlag
        
    # loop to find root
    while not done:
            iter_num += 1
            xr = (lb + ub) / 2
            fun_xr = sb.F_net(np.array([xr]), np.array([x]), sb.rover, sb.planet, crr)
            
            # sub in xr
            if fun_l * fun_xr < 0:
                 ub = xr
            else:
                 lb = xr

            err_est = abs((lb - ub) / (lb + ub))

            # check limits
            if err_est < err_max:
                 done = True
                 root = xr
                 exitFlag = 1 # all good

            if iter_num >= iter_max:
                 done = True
                 root = xr
                 exitFlag = 2 # iter_max reached
    
    # find velocity with gear_ratio
    speed_reducer = sb.rover['wheel_assembly']['speed_reducer']
    gear_ratio = sb.get_gear_ratio(speed_reducer)
    velocity = root * sb.rover['wheel_assembly']['wheel']['radius'] / gear_ratio
    
    return velocity, err_est, iter_num, exitFlag



for angle in range(len(slope_array_deg)):
    velocity, err_est, iter_num, exitFlag = bisection(slope_array_deg[angle])
    v_max.append(velocity)
    iter_num_list.append(iter_num)
    err_est_list.append(err_est)
    exitFlag_list.append(exitFlag)

plt.plot(slope_array_deg, v_max, color='r')
plt.title('Max rover speed vs Terrain slope')
plt.ylabel('Maximum Velocity (m/s)')
plt.xlabel('Slope of Terrain (deg)')
plt.show()
