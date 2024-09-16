#By submitting this assignment, I agree to the following:
#   "Aggies do not lie, cheat, or steal, or tolerate those who do."
#   "I have not given or received any unauthorized aid on this assignment."
#
# Name:        Jackson Walton
# Section:      505
# Assignment:   THE ASSIGNMENT NUMBER (e.g. Lab 1b-2)
# Date:         DAY MONTH 2024

import matplotlib.pyplot as plt
from subfunctions import *
import numpy as np

Crr_array = np.linspace(0.01, 0.4, 25)
v_max = []

def bisection(Crr):
    lb = 0
    ub = 3.8
    err_max = 1e-6
    iter_max = 50
    
    done = False
    numIter = 0
    
    root = np.nan      
    err = np.nan
    
    fl = F_net(np.array([lb]),np.array([0]),rover, planet, Crr)
    fu = F_net(np.array([ub]),np.array([0]),rover, planet, Crr)
    
    if fl * fu > 0:    
        done = True
    elif fl * fu == 0:  
        err = 0
        done = True
        
        if abs(fl) <= abs(fu):  
            root = lb
        else:
            root = ub
        
    while not done:
        
        numIter += 1
        
        xr = (ub + lb) / 2
        fr = F_net(np.array([xr]),np.array([0]),rover, planet, Crr)
        
        if fr is np.nan:
            done = True
            
        if fl * fr == 0:    
            done = True
            err = 0
            root = xr
        elif fl * fr > 0:   
            lb = xr
        else:               
            ub = xr
            
        err = 100 * abs((ub - lb) / xr)     
        
        if err <= err_max:      
            root = xr
            done = True
            break
                  
        if numIter >= iter_max:     
            root = xr
            done = True
            break  
    
    speed_reducer = rover['wheel_assembly']['speed_reducer']
    gear_ratio = get_gear_ratio(speed_reducer)
    root = root * rover['wheel_assembly']['wheel']['radius'] / gear_ratio 
    
    return root


for i in range(len(Crr_array)):
    v_i = bisection(Crr_array[i])
    v_max.append(v_i)
    
    
plt.plot(Crr_array, v_max)
plt.title("Max Rover Speed vs Coefficient of Rolling Resistance")
plt.ylabel('Max Velocity (m/s)')   
plt.xlabel('Crr Array (unitless)') 

