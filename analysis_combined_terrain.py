#By submitting this assignment, I agree to the following:
#   "Aggies do not lie, cheat, or steal, or tolerate those who do."
#   "I have not given or received any unauthorized aid on this assignment."
#
# Name:        Jackson Walton
# Section:      505
# Assignment:   THE ASSIGNMENT NUMBER (e.g. Lab 1b-2)
# Date:         DAY MONTH 2024

import matplotlib.pyplot as plt
import subfunctions as sb
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


Crr_array = np.linspace(0.01, 0.4, 25)
slope_array_deg = np.linspace(-10, 35, 25)
CRR, SLOPE = np.meshgrid(Crr_array, slope_array_deg)
VMAX = np.zeros(np.shape(CRR), dtype = float)
N = np.shape(CRR)[0]

def bisection(slope_array_deg, Crr): #find roots for F_net (F=0, v=max)

    lb = 0
    ub = 3.8 # m/s
    err_max = 1e-6
    iter_max = 50

    # y-values for lower and upper bound
    fun_l = sb.F_net(np.array([lb]), np.array([slope_array_deg]), sb.rover, sb.planet, Crr)
    fun_u = sb.F_net(np.array([ub]), np.array([slope_array_deg]), sb.rover, sb.planet, Crr)

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
            fun_xr = sb.F_net(np.array([xr]), np.array([slope_array_deg]), sb.rover, sb.planet, Crr)
            
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
    
    return root

# def bisection(slope_array_deg, Crr, rover, planet):
#     lb = 0
#     ub = 3.8
#     err_max = 1e-6
#     iter_max = 30
    
#     done = False
#     numIter = 0
    
#     root = np.nan      
#     app_err = np.nan
    
#     fl = F_net(np.array([lb]), np.array([slope_array_deg]), rover, planet, Crr)
#     fu = F_net(np.array([ub]), np.array([slope_array_deg]), rover, planet, Crr)
    
#     if fl * fu > 0:     
#         done = True
#     elif fl * fu == 0: 
#         done = True
#         app_err = 0.0
        
#         if abs(fl) <= abs(fu):  
#             root = lb
#         else:
#             root = ub

#     while not done:
#         numIter += 1
#         xr = (ub + lb) / 2
#         fc = F_net(np.array([xr]), np.array([slope_array_deg]), rover, planet, Crr)
        
#         if fc is np.nan:
#             done = True
#             break
        
#         if fl * fc == 0.0:    
#             done = True
#             app_err = 0.0
#             root = xr
#             break   
#         elif fl * fc > 0:   
#             lb = xr
#         else:               
#             ub = xr
            
#         app_err = abs((ub - lb) / xr)     
        
#         if app_err <= err_max:      
#             done = True
#             root = xr
#             break       
            
#         if numIter >= iter_max:     
#             done = True
#             root = xr   
#             break       
    
    
#     speed_reducer = rover['wheel_assembly']['speed_reducer']
#     gear_ratio = get_gear_ratio(speed_reducer)
#     root = root * rover['wheel_assembly']['wheel']['radius'] / gear_ratio
    
#     return root


for i in range(N):
    for j in range(N):
        Crr_sample = float(CRR[i,j])
        slope_sample = float(SLOPE[i,j])
        VMAX[i,j] = bisection(slope_sample, Crr_sample)


figure = plt.figure()
ax = Axes3D(figure, elev = 15, azim = 30) #where N1 and N2 will control the 3D view
figure.add_axes(ax)
ax.plot_surface(CRR, SLOPE, VMAX)
ax.set_xlabel("Crr_array (unitless)")
ax.set_ylabel("slope_array_deg (deg)")
ax.set_zlabel("v_max (m/s)")
ax.set_title("Speed of Rover vs Rolling Resistance and Terrain Slope")
ax.show()