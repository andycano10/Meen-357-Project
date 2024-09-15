# -*- coding: utf-8 -*-
#By submitting this assignment, I agree to the following:
#   "Aggies do not lie, cheat, or steal, or tolerate those who do."
#   "I have not given or received any unauthorized aid on this assignment."
#
# Name:        Jackson Walton, Andy Cano-Avila, Catalina Fossati
# Section:      505
# Assignment:   THE ASSIGNMENT NUMBER (e.g. Lab 1b-2)
# Date:         DAY MONTH 2024

import math
import numpy as np

######################## Dictionaries ########################

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
    'mass': 5.0}

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
    'power_subsys': power_subsys}


######################## Functions ########################

def get_mass(rover):

    # Check that the input is a dict
    if type(rover) != dict:
        raise Exception('Argument must be a dict')
# initialize
    m = 0
    m += (rover['wheel_assembly']['wheel']['mass']) * 6
    m += (rover['wheel_assembly']['motor']['mass']) * 6
    m += (rover['wheel_assembly']['speed_reducer']['mass']) * 6
    m += rover['chassis']['mass']
    m += rover['science_payload']['mass']
    m += rover['power_subsys']['mass']

    return m



def get_gear_ratio(speed_reducer):

    # Check that the input is a dict
    if type(speed_reducer) != dict:
        raise Exception('Argument must be a dict')

    # Check the type of the dict, not case sensitive
    if speed_reducer['type'].lower() != 'reverted':
        raise Exception('The type of speed reducer is invalid.')

    # Main function
    diam1 = speed_reducer['diam_pinion']
    diam2 = speed_reducer['diam_gear']

    Ng = (diam2/diam1)**2

    return Ng





def tau_dcmotor(omega, motor):

    # Exceptions
    if not isinstance(omega, np.ndarray):
        raise Exception('Please input a scalar or vector. No matrices.')

    if type(motor) != dict:
        raise Exception('Motor input must be a dict')

    for i in omega:
        try:
            n = float(i)
        except:
            raise Exception('Please input a scalar or vector. No matrices.')

    T_s = motor['torque_stall']
    T_n = motor['torque_noload']
    w_n = motor['speed_noload']
    w = omega
    tau = []
   
    for i in w:
        if i > w_n:
            tau += [0]
        elif i < 0:
            tau += [T_s]
        else:
            tau += [T_s - ((T_s - T_n)/(w_n)) * i]
    return tau

def F_gravity(terrain_angle, rover, planet):

  # Exceptions
  if not isinstance(terrain_angle, np.ndarray) and not np.isscalar(terrain_angle):
      raise Exception('Terrain angle input must be a vector or a scalar')

  if not (-75 <= np.min(terrain_angle) <= 75) or not (-75 <= np.max(terrain_angle) <= 75):
      raise Exception('Inputs for terrain angles must be between -75 and +75 degrees')

  if not isinstance(rover, dict) or not isinstance(planet, dict):
      raise Exception('Inputs for rover and planet must be dicts')

  rov_mass = get_mass(rover)

  ga = planet['g']

  Fgt = -1*rov_mass * ga * np.sin(np.radians(terrain_angle))

  return Fgt



def F_rolling(omega, terrain_angle, rover, planet, Crr):

    for i in omega:
        try:
            n = float(i)
        except:
            raise Exception('Please input a scalar or vector. No matrices.')
           
    for i in terrain_angle:
        try:
            n = float(i)
        except:
            raise Exception('Please input a scalar or vector. No matrices.')
           
    if len(omega) != len(terrain_angle):
        raise Exception("Both vectors should be equivalent sizes")
    if not (-75 <= np.min(terrain_angle) <= 75) or not (-75 <= np.max(terrain_angle) <= 75):
        raise Exception('Inputs for terrain angles must be between -75 and +75 degrees')
    if Crr <= 0:
        raise Exception("The coefficient of rolling resistance must be positive")
    if (type(rover) != dict) or (type(planet) != dict):
        raise Exception('Inputs for rover and planet must be dicts')

    rov_mass = get_mass(rover)
    ga = planet['g']
    Frr = []
    j = 0
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
   
    for i in terrain_angle:
        w = omega[j] / Ng
        F_n = rov_mass * ga * np.cos(np.radians(i)) * Crr
        v_rov = rover['wheel_assembly']['wheel']['radius'] * w
        Frr.append(math.erf(40 * v_rov) * -1*F_n)
        j += 1

    return Frr



def F_drive(omega, rover):

  if not isinstance(omega, np.ndarray) and not np.isscalar(omega):
      raise Exception('The shaft speed must be a vector or scalar')

  if not isinstance(rover, dict):
      raise Exception('Rover input type must be a dict')


  mot = rover['wheel_assembly']['motor']
  Fd = []
  T_in = tau_dcmotor(omega, mot)

  Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])

  T_out = T_in

  wheel_rad = rover['wheel_assembly']['wheel']['radius']
  for i in T_out:
      Fd += [(6 * i * Ng)/ wheel_rad]
  return Fd



def F_net(omega, terrain_angle, rover, planet, Crr):

    # Exceptions
    if not np.isscalar(omega) and not isinstance(omega, np.ndarray):
        raise Exception('Motor shaft speed must be a scalar or a vector')

    if not np.isscalar(terrain_angle) and not isinstance(terrain_angle, np.ndarray):
        raise Exception('The terrain angle must be a scalar or a vector')

    for i in terrain_angle:
        if i < -75 or i > 75:
            raise Exception('All terrain angles must be between -75 and +75 degrees')

    if len(omega) != len(terrain_angle):
        raise Exception("Omega and the terrain angle should be equivalent lengths")

    if not isinstance(rover, dict) or not isinstance(planet, dict):
        raise Exception('Rover input must be a dict')

    if not np.isscalar(Crr) or Crr <= 0:
        raise Exception('This value must be a positive scalar')

    Fg = F_gravity(terrain_angle, rover, planet)
    Fd = F_drive(omega, rover)
    Fr = F_rolling(omega, terrain_angle, rover, planet, Crr)
    Fn = Fg + Fd + Fr

    return Fn
