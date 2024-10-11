
import matplotlib.pyplot as plt
import subfunctions as sb
import numpy as np

Crr_array = np.linspace(0.01, 0.4, 25)
v_max = []

def bisection(Crr):

    lb = 0
    ub = 3.8 # m/s
    err_max = 1e-6
    iter_max = 50

    # y-values for lower and upper bound
    fun_l = sb.F_net(np.array([lb]), np.array([0]), sb.rover, sb.planet, Crr)
    fun_u = sb.F_net(np.array([ub]), np.array([0]), sb.rover, sb.planet, Crr)

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
    
    elif fun_l * fun_u == 0: 
        done = True  
        
    # loop to find root
    while not done:
            iter_num += 1
            xr = (lb + ub) / 2
            fun_xr = sb.F_net(np.array([xr]), np.array([0]), sb.rover, sb.planet, Crr)
            
            if fun_l * fun_xr == 0.0:    
                done = True

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
    root = root * sb.rover['wheel_assembly']['wheel']['radius'] / gear_ratio
    
    return root, err_est, iter_num, exitFlag

for i in range(len(Crr_array)):
    v_i, err_est, iter_num, exitFlag = bisection(Crr_array[i])
    v_max.append(v_i)
    
    
plt.plot(Crr_array, v_max)
plt.title("Max Rover Speed vs Coefficient of Rolling Resistance")
plt.ylabel('Max Velocity (m/s)')   
plt.xlabel('Crr Array (unitless)') 