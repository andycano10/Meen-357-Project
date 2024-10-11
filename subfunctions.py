"""###########################################################################
#   This file contains subfunctions for Phase 1 of the TAMU MEEN 357 project
#
#   Created by: MEEN 357 Instructional Team
#   Last Modified: 4 October 2023
###########################################################################"""

import math
import numpy as np
from statistics import mean
from scipy.interpolate import *
from scipy.integrate import *

wheel = {
    'radius': 0.3,
    'mass': 1.0 }

speed_reducer = {
    'type': 'reverted',
    'diam_pinion': 0.04,
    'diam_gear': 0.07,
    'mass': 1.5}

motor = {
    'torque_stall': 170,
    'torque_noload': 0,
    'speed_noload': 3.80,
    'mass': 5.0,
    'effcy_tau':np.array([0, 10, 20, 40, 70, 165]),
    'effcy': np.array([0, 0.55, 0.75, 0.71, 0.50, 0.05])}

chassis = {'mass': 659}
science_payload = {'mass': 75}
power_subsys = {'mass': 90}

planet = {'g': 3.72}

wheel_assembly = {
    'wheel': wheel,
    'speed_reducer': speed_reducer,
    'motor': motor}

rover = {
    'wheel_assembly': wheel_assembly,
    'chassis': chassis,
    'science_payload': science_payload,
    'power_subsys': power_subsys,
    'telemetry': None}

def get_mass(rover):
    """
    Inputs:  rover:  dict      Data structure containing rover parameters
    
    Outputs:     m:  scalar    Rover mass [kg].
    """
    
    # Check that the input is a dict
    if type(rover) != dict:
        raise Exception('Input must be a dict')
    
    # add up mass of chassis, power subsystem, science payload, 
    # and components from all six wheel assemblies
    m = rover['chassis']['mass'] \
        + rover['power_subsys']['mass'] \
        + rover['science_payload']['mass'] \
        + 6*rover['wheel_assembly']['motor']['mass'] \
        + 6*rover['wheel_assembly']['speed_reducer']['mass'] \
        + 6*rover['wheel_assembly']['wheel']['mass'] \
    
    return m


def get_gear_ratio(speed_reducer):
    """
    Inputs:  speed_reducer:  dict      Data dictionary specifying speed
                                        reducer parameters
    Outputs:            Ng:  scalar    Speed ratio from input pinion shaft
                                        to output gear shaft. Unitless.
    """
    
    # Check that the input is a dict
    if type(speed_reducer) != dict:
        raise Exception('Input must be a dict')
    
    # Check 'type' field (not case sensitive)
    if speed_reducer['type'].lower() != 'reverted':
        raise Exception('The speed reducer type is not recognized.')
    
    # Main code
    d1 = speed_reducer['diam_pinion']
    d2 = speed_reducer['diam_gear']
    
    Ng = (d2/d1)**2
    
    return Ng


def tau_dcmotor(omega, motor):
    """
    Inputs:  omega:  numpy array      Motor shaft speed [rad/s]
             motor:  dict             Data dictionary specifying motor parameters
    Outputs:   tau:  numpy array      Torque at motor shaft [Nm].  Return argument
                                      is same size as first input argument.
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(motor) != dict:
        raise Exception('Second input must be a dict')
        
    # Main code
    tau_s    = motor['torque_stall']
    tau_nl   = motor['torque_noload']
    omega_nl = motor['speed_noload']
    
    # initialize
    tau = np.zeros(len(omega),dtype = float)
    for ii in range(len(omega)):
        if omega[ii] >= 0 and omega[ii] <= omega_nl:
            tau[ii] = tau_s - (tau_s-tau_nl)/omega_nl *omega[ii]
        elif omega[ii] < 0:
            tau[ii] = tau_s
        elif omega[ii] > omega_nl:
            tau[ii] = 0
        
    return tau
    
    


def F_rolling(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  numpy array     Motor shaft speed [rad/s]
              terrain_angle:  numpy array     Array of terrain angles [deg]
                      rover:  dict            Data structure specifying rover 
                                              parameters
                    planet:  dict            Data dictionary specifying planetary 
                                              parameters
                        Crr:  scalar          Value of rolling resistance coefficient
                                              [-]
    
    Outputs:           Frr:  numpy array     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    r = rover['wheel_assembly']['wheel']['radius']
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    v_rover = r*omega/Ng
    
    Fn = np.array([m*g*math.cos(math.radians(x)) for x in terrain_angle],dtype=float) # normal force
    Frr_simple = -Crr*Fn # simple rolling resistance
    
    Frr = np.array([math.erf(40*v_rover[ii]) * Frr_simple[ii] for ii in range(len(v_rover))], dtype = float)
    
    return Frr


def F_gravity(terrain_angle, rover, planet):
    """
    Inputs:  terrain_angle:  numpy array   Array of terrain angles [deg]
                     rover:  dict          Data structure specifying rover 
                                            parameters
                    planet:  dict          Data dictionary specifying planetary 
                                            parameters
    
    Outputs:           Fgt:  numpy array   Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that values of the first input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the first input must be between -75 degrees and +75 degrees')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Check that the third input is a dict
    if type(planet) != dict:
        raise Exception('Third input must be a dict')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    
    Fgt = np.array([-m*g*math.sin(math.radians(x)) for x in terrain_angle], dtype = float)
        
    return Fgt


def F_drive(omega, rover):
    """
    Inputs:  omega:  numpy array   Array of motor shaft speeds [rad/s]
             rover:  dict          Data dictionary specifying rover parameters
    
    Outputs:    Fd:  numpy array   Array of drive forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Main code
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor'])
    tau_out = tau*Ng
    
    r = rover['wheel_assembly']['wheel']['radius']
    
    # Drive force for one wheel
    Fd_wheel = tau_out/r 
    
    # Drive force for all six wheels
    Fd = 6*Fd_wheel
    
    return Fd


def F_net(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  list     Motor shaft speed [rad/s]
              terrain_angle:  list     Array of terrain angles [deg]
                      rover:  dict     Data structure specifying rover 
                                      parameters
                     planet:  dict     Data dictionary specifying planetary 
                                      parameters
                        Crr:  scalar   Value of rolling resistance coefficient
                                      [-]
    
    Outputs:           Fnet:  list     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
    # if (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
    
    # Main Code
    Fd = F_drive(omega, rover)
    Frr = F_rolling(omega, terrain_angle, rover, planet, Crr)
    Fg = F_gravity(terrain_angle, rover, planet)
    
    Fnet = Fd + Frr + Fg # signs are handled in individual functions
    
    return Fnet


# #%% Hint for students
# omega = 1 # rad/s (motor shaft speed)
# angle = 5 # degrees (terrain angle)
# Crr = 0.1

# from define_rovers import define_rover_1
# rover, planet = define_rover_1() # this is what is given in Appendix A and B

# Fd = F_drive(omega, rover)
# Frr = F_rolling(omega, angle, rover, planet, Crr)
# Fg = F_gravity(angle, rover, planet)
# Fnet = F_net(omega, angle, rover, planet, Crr)


################################### Phase 2 ###################################
def experiment1():
    
    experiment = {'time_range' : np.array([0,20000]),
                  'initial_conditions' : np.array([0.3025,0]),
                  'alpha_dist' : np.array([0, 100, 200, 300, 400, 500, 600, \
                                           700, 800, 900, 1000]),
                  'alpha_deg' : np.array([11.509, 2.032, 7.182, 2.478, \
                                        5.511, 10.981, 5.601, -0.184, \
                                        0.714, 4.151, 4.042]),
                  'Crr' : 0.1}
    
    
    # Below are default values for example only:
    end_event = {'max_distance' : 50,
                 'max_time' : 5000,
                 'min_velocity' : 0.01}
    
    return experiment, end_event




def end_of_mission_event(end_event):
    """
    Defines an event that terminates the mission simulation. Mission is over
    when rover reaches a certain distance, has moved for a maximum simulation 
    time or has reached a minimum velocity.            
    """
    
    mission_distance = end_event['max_distance']
    mission_max_time = end_event['max_time']
    mission_min_velocity = end_event['min_velocity']
    
    # Assume that y[1] is the distance traveled
    distance_left = lambda t,y: mission_distance - y[1]
    distance_left.terminal = True
    
    time_left = lambda t,y: mission_max_time - t
    time_left.terminal = True
    
    velocity_threshold = lambda t,y: y[0] - mission_min_velocity;
    velocity_threshold.terminal = True
    velocity_threshold.direction = -1
    
    # terminal indicates whether any of the conditions can lead to the
    # termination of the ODE solver. In this case all conditions can terminate
    # the simulation independently.
    
    # direction indicates whether the direction along which the different
    # conditions is reached matters or does not matter. In this case, only
    # the direction in which the velocity treshold is arrived at matters
    # (negative)
    
    events = [distance_left, time_left, velocity_threshold]
    
    return events

def motorW(v, rover):
   
    '''
    
    This function is used to calculate the rotational speed of the 
    motor shaft (in rad/s) using the rover dictionary. and the translational 
    velocity of the rover. 
    
    Inputs: 
        v : a scalar or a numpy array that represents translational velocity
        rover : a rover dictionary
    
    Outputs:
        w : the rotational speed of the motor shaft (in rad/s)
        
    '''  
    #Check errors
    if not isinstance(v, np.ndarray) and not np.isscalar(v):
        raise Exception('Your input for velocity must be an array or scalar.')
        
    if not isinstance(rover, dict):
        raise Exception('Th input for rover must be a dict.')
    # Gear ratio function
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    # Establish radius variable
    rad = rover['wheel_assembly']['wheel']['radius']
    
    # Calc rotational speed for scalars
    if np.isscalar(v):
        w = (v / rad) * Ng
        
    # Calc rotational speeds of each index and add them to final array
    if isinstance(v, np.ndarray):
        w = []
        for i in v:
            w += [(i / rad) * Ng]
            
    return w

def rover_dynamics(t, y, rover, planet, experiment):
    
    '''
    
    This function is used to determine the derivative of the state vector for 
    the rover given current conditions.
    
    Inputs:
        t : a scalar representing a time sample
        y : a 1D numpy array that represents the state vector
        rover : a rover dictionary
        planet : a planet dictionary
        experiment : an experiment dictionary
        
    Outputs:
        dydt : derivative of the state vector
    
    '''

# Check exceptions    
        
    if not isinstance(y, np.ndarray):
        raise Exception('The input for y must be a 1D array')
        
    if not np.isscalar(t):
        raise Exception('The input for t must be scalar.')
        
    if not isinstance(rover, dict):
        raise Exception('The input for rover must be a dict')
        
    if not isinstance(planet, dict):
        raise Exception('The input for planet must be a dict')
    
    if not isinstance(experiment, dict):
        raise Exception('The input for experiment must be a dict')
        
    # Call mass function
    mass = get_mass(rover)
    
    # Create rotational speed array
    w = np.array([motorW(y[0], rover)])
    
    # Assign given values
    Crr = experiment['Crr']   
    alpha_deg = experiment['alpha_deg']
    alpha_dist = experiment['alpha_dist']

    # Interpolate terrain angle at different angles and distances
    alpha_fun = interp1d(alpha_dist, alpha_deg, kind = 'cubic', fill_value='extrapolate')
    
    terr_angle = np.array([alpha_fun(y[1])])
    
    Fnet = float(F_net(w, terr_angle, rover, planet, Crr))
    # This calculates the acceleration using the net force and mass
    acc = Fnet / mass
    
    dydt = np.array([acc, y[0]])
    
    return dydt


def mechpower(v, rover):

    '''
    
    This function is used to calculate the power output at an instantaneous 
    moment by a singlular DC motor at each point within a known velocity profile.
    
    Inputs: 
        v : a scalar or numpy array that represents translational velocity
        rover : a rover dictionary
    
    Outputs:
        w : the rotational speed of the motor shaft (in rad/s)
        
    '''
    
    # Check exceptions
    if not isinstance(v, np.ndarray) and not np.isscalar(v):
        raise Exception('Input for velocity must be an array or scalar')
                    
    if not isinstance(rover, dict):
        raise Exception('Input for rover must be a dict')
        
    mot = rover['wheel_assembly']['motor']
    
    # Calc power using rotational speed and torque
    if isinstance(v, np.ndarray):
        w = []
        for i in v:
            w += [motorW(i, rover)]
            
        w_arr = np.array(w)        
        tau = tau_dcmotor(w_arr, mot)       
        tau_arr = np.array(tau)       
        P = tau_arr * w_arr
        
    # Calc power for scalars
    if np.isscalar(v):
        
        w = np.array([motorW(v, rover)])     
        tau = np.array([tau_dcmotor(w, mot)])    
        P = float(tau * w)
    
    return P

def battenergy(t, v, rover):
   
    '''
    
    This function is used to calculate the total electrical energy consumed 
    from the rover battery pack over a simulation profile, defined as 
    time-velocity pairs. It assumes all six motors are driven from 
    the same battery pack.
    
    Inputs:
        t : a 1D numpy array that represents samples of time from a rover simulation
        v : a 1D numpy array that represents velocity data from a rover simulation
        rover : a rover dictionary
        
    Outputs:
        E : a scalar representing the total electrical energy consumed throughout the simulation profile
    
    ''' 
   
    if not isinstance(v, np.ndarray):
        raise Exception('Input for velocity must be a 1D array')
        
    if not isinstance(t, np.ndarray):
        raise Exception('Input for time must be a 1D array')     
        
    if len(t) != len(v):
        raise Exception('The time and velocity inputs must be the same length')
        
    if not isinstance(rover, dict):
        raise Exception('The inpit for rover must be a dict')
 
    w = np.array(motorW(v, rover))
    P = mechpower(v, rover)
    mot = rover['wheel_assembly']['motor']  
    tau = tau_dcmotor(w, mot)
    
    effcy_tau = rover['wheel_assembly']['motor']['effcy_tau']
    effcy = rover['wheel_assembly']['motor']['effcy']
    
    # Interpolate efficiency
    effcy_fun = interp1d(effcy_tau, effcy, kind = 'cubic') 
    effcy_val = effcy_fun(tau)
    
    # Calc battery power with efficiency and mech power
    p_battery = P / effcy_val
    # Integrate power w/ respect to time
    E_single = np.trapz(p_battery, t)
    
    E = 6 * E_single
    
    return E

def simulate_rover(rover, planet, experiment, end_event):
    
    '''
    
    This function is used to integrate the trajectory of the rover according 
    to terrain and initial conditions and uses this data to populate the 
    'telemetry' section within the rover dictionary.
    
    Inputs:
        rover : a rover dictionary
        planet : a planet dictionary
        experiment : an experiment dictionary 
        end_event : an end_event dictionary
        
    Outputs:
        rover : an updated rover dictionary with the 'telemetry' data
    
    '''
    # Check exceptions
    if not isinstance(rover, dict):
        raise Exception('The input for rover must be a dict')
    
    if not isinstance(planet, dict):
        raise Exception('The input for planet must be a dict')
    
    if not isinstance(experiment, dict):
        raise Exception('The input for experiment must be a dict')
    
    if not isinstance(end_event, dict):
        raise Exception('The input for end_event must be a dict')
    
    
    # Create rover dynamics function
    dyn_func = lambda t,y: rover_dynamics(t, y, rover, planet, experiment)     
    timespan = experiment['time_range']    
    y_0 = experiment['initial_conditions'].ravel() 
    
    event = end_of_mission_event(end_event) 
    # Solve IVP
    solut = solve_ivp(dyn_func, timespan, y_0, method = 'BDF', event=event) 
    
    
    time = solut.t    
    completion_t = solut.t[-1] 
    vel = solut.y[0,:] 
    pos = solut.y[1,:]  
    dist_trav = solut.y[1,-1]  
    max_vel = max(solut.y[0,:])    
    avg_vel = mean(solut.y[0,:])   
    power = mechpower(solut.y[0,:], rover)   
    bat_energy = battenergy(solut.t,solut.y[0,:],rover)  
    energy_per_dist = battenergy(solut.t,solut.y[0,:],rover)/solut.y[1,-1]

    # Update telemetry key to include new data
    telemetry = {'Time':time,
                 'completion_time':completion_t,
                 'velocity':vel,
                 'position':pos,
                 'distance_traveled':dist_trav,
                 'max_velocity':max_vel,
                 'average_velocity':avg_vel,
                 'power':power,
                 'battery_energy':bat_energy,
                 'energy_per_distance':energy_per_dist}
    
    # Add telemetry info to rover dict
    rover['telemetry'] = telemetry
      
    return rover
