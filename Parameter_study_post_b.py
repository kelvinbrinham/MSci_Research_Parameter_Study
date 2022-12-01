#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:54:09 2022

@author: Kelvin Brinham
"""

'''
This parameter study searches for viscous instability over a range of mass accretion rates and radii for steady circumstellar discs.

It consists of 3 scripts, this is part b

Part b finds the temperature solutions, snowline radii, optical depth solutions, surface density solutions and alpha solutions for 
a range of mass accretion rates.

The data produced here is saved, this allows it to be used in script c without re-running this potentially long calculation.

NB: This parameter study uses cgs units unless stated otherwise
NB: Interior/Exterior means at a lower/higher radius than the snowline
'''

#%%
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize

import Parameter_study_post_a



'''
Get Temperature, Sigma and Snowline radii data for multiple Mdot
'''

Radii_inner_index = [] #Indices of snowline inner edge radii in the radius array for each Mdot
Radii_outer_index = [] #Indices of snowline outer edge radii in the radius array for each Mdot

Radii_inner = [] #Inner radius of snowline for each Mdot
Radii_middle = [] #Middle radius of snowline for each Mdot
Radii_outer = [] #Outer radius of snowline for each Mdot

Radii_snowline = [Radii_inner, Radii_outer] # Defined for ease


'''
Useful key to keep track of data in arrays

Temperature Key:
T[0] = [np.array([T_soln1]), np.array([T_soln2]), np.array([T_soln3])] for Mdot1
T[1] = [np.array([T_soln1]), np.array([T_soln2]), np.array([T_soln3])] for Mdot2
T[0][0] = np.array([T_soln1]) for Mdot1

Sigma Key:
Sigma_data[0] = [np.array([Sigma_soln1]), np.array([Sigma_soln2]), np.array([Sigma_soln3])] for Mdot1
Sigma_data[0][0] = np.array([Sigma_soln1]) for Mdot1

Tau Key:
Tau_data[0] = [[[Tau_soln1]], [[Tau_soln2]], [[Tau_soln3]]]  for Mdot1
Tau_data[0][0] = [[Tau_soln1]] for Mdot1
'''

T = [] # Temperature solutions (see above key) e.g [[[T_Mdot1_soln1], [T_Mdot1_soln2], [T_Mdot1_soln3]], ..., [[T_Mdotn_soln1], [T_Mdotn_soln2], [T_Mdotn_soln3]]]

Sigma_data = [] #Surface density solutions(see above key)
Tau_data = [] #Optical depth solutions(see above key)
alpha_data = [] #Alpha solutions(see above key)

t = 0 # Keeps track of which Mdot solver is on, useful for knowing how the solver is progressing

for Mdot_current in Mdot_array:
    
    t += 1
    print(t)
    
    # Calculate Temperature 
    Temperature_data = np.array(T_Mdot_combination(R, Mdot_current)) # np.array([[], [], []]) Array of 3 lists, one for each Temp. solution all for 1 Mdot
    T.append(Temperature_data)
    Temperature_data_inner = abs(Temperature_data[0] - T_snowline_start) # Use first solution for inner temperature of snowline (because jumps down at start)
    Temperature_data_outer = abs(Temperature_data[2] - T_snowline_end) # Use third solution for outer temperature of snowline (because jumps down at end)


    #Calculate Snowline Radii
    R_inner_index = np.argmin(Temperature_data_inner)
    R_outer_index = np.argmin(Temperature_data_outer)

    
    R_inner = R[R_inner_index]
    R_outer = R[R_outer_index]
    R_middle = 0.5*(R_outer - R_inner) + R_inner
    # NB: IF no snowline present in array of R then R_inner and R_outer return R[0] or R[-1]
    Radii_inner_index.append(R_inner_index)
    Radii_outer_index.append(R_outer_index)
    
    Radii_inner.append(R_inner)
    Radii_middle.append(R_middle)
    Radii_outer.append(R_outer)
    
    #Caclulate Tau
    Tau_data_current = tau_data(R)
    Tau_data.append(Tau_data_current)
    
    
    #Calculate Sigma and alpha_gt
    sigma_data_single_Mdot = [] #[[], [], []]
    alpha_gt_data_single_Mdot = [] #[[], [], []]
    
    for i in range(3):
        
        sigma = Sigma(Temperature_data[i], R)
        
        sigma_data_single_Mdot.append(sigma)
        
        
        alpha = alpha_gt_data(Temperature_data[i], sigma, Tau_data_current[i][0], R)
        
        alpha_gt_data_single_Mdot.append(alpha)
        
        
    Sigma_data.append(sigma_data_single_Mdot)
    alpha_data.append(alpha_gt_data_single_Mdot)
    
    

Radius_array = R[:(Radii_outer_index[0] + 1)][Radii_inner_index[-1]:] #All the radii within snowlines 
# (i.e. All the radii from the inner edge of the innermost snowline to the outer edge of the outermost snowline)

        

#%%
'''
Saving the above data 
'''
#Transform to numpy arrays
G = np.asarray(G)
M = np.asarray(M)
Mdot_array = np.asarray(Mdot_array)
Q0 = np.asarray(Q0)
R = np.asarray(R)
Radii_inner = np.asarray(Radii_inner)
Radii_inner_index = np.asarray(Radii_inner_index)
Radii_middle = np.asarray(Radii_middle)
Radii_outer = np.asarray(Radii_outer)
Radii_outer_index = np.asarray(Radii_outer_index)
Radii_snowline = np.asarray(Radii_snowline)
Radius_array = np.asarray(Radius_array)
Sigma_data = np.asarray(Sigma_data)
T = np.asarray(T)
T_snowline_end = np.asarray(T_snowline_end)
Tau_data = np.asarray(Tau_data)
Temperature_data = np.asarray(Temperature_data)
Temperature_data_inner = np.asarray(Temperature_data_inner)
Temperature_data_outer = np.asarray(Temperature_data_outer)
Tirr = np.asarray(Tirr)
alpha_data = np.asarray(alpha_data)
alpha_frag = np.asarray(alpha_frag)
au = np.asarray(au)
beta = np.asarray(beta)
dT = np.asarray(dT)
kappa0 = np.asarray(kappa0)
kappa_drop = np.asarray(kappa_drop)
kb = np.asarray(kb)
mu = np.asarray(mu)
sb = np.asarray(sb)
sm = np.asarray(sm)
gamma = np.asarray(gamma)
Tirr = np.asarray(Tirr)

#Saving data
np.save('G', G)
np.save('M', M)
np.save('Mdot_array', Mdot_array)
np.save('Q0', Q0)
np.save('R', R)
np.save('Radii_inner', Radii_inner)
np.save('Radii_inner_index', Radii_inner_index)
np.save('Radii_middle', Radii_middle)
np.save('Radii_outer', Radii_outer)
np.save('Radii_outer_index', Radii_outer_index)
np.save('Radii_snowline', Radii_snowline)
np.save('Radius_array', Radius_array)
np.save('Sigma_data', Sigma_data)
np.save('T', T)
np.save('T_snowline_end', T_snowline_end)
np.save('Tau_data', Tau_data)
np.save('Temperature_data', Temperature_data)
np.save('Temperature_data_inner', Temperature_data_inner)
np.save('Temperature_data_outer', Temperature_data_outer)
np.save('Tirr', Tirr)
np.save('alpha_data', alpha_data)
np.save('alpha_frag', alpha_frag)
np.save('au', au)
np.save('beta', beta)
np.save('dT', dT)
np.save('kappa0', kappa0)
np.save('kappa_drop', kappa_drop)
np.save('kb', kb)
np.save('mu', mu)
np.save('sb', sb)
np.save('sm', sm)
np.save('gamma', gamma)
np.save('Tirr', Tirr)

