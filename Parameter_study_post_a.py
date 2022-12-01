#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:15:09 2022

@author: Kelvin Brinham
"""

'''
This parameter study searches for viscous instability over a range of mass accretion rates and radii for steady circumstellar discs.

It consists of 3 scripts, this is part a

Part a sets the constants and defines the functions needed

NB: This parameter study uses cgs units unless stated otherwise
NB: Interior/Exterior means at a lower/higher radius than the snowline
'''

#%%

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize

#%%
# Definining constants       

gamma = 5/3 #Ratio of specific heats
G = 6.67e-8 #Gravitational constant
kb = 1.38e-16 #Boltzmann constant
sb = 5.6704e-5 #Stefan-Boltzmann constant
au = 1.495978707e13 #1 au in cm
sm = (2e33)/(365*24*60**2) #1 Solar mass per year in g/s


#Parameters

mu = 2.35*1.67e-24 #Mean particle mass in g, doesnt include solids, 2.35 from H2 and He ratios and small amounts
#mu = 3.38*1.67e-24 #Rafikov mu
#mu = 1.18*1.67e-24 #Rafikov value from fidicual values


beta = 2.0 

kappa0 = 1.41e-5 #Minimum possible for CO ice

Q0 = 1.5 #Critical Toomre factor 

M = 2e33 # Stellar mass


alpha_frag = 0.1 #Fragmentation threshold alpha

Tirr = 12 # Irradiation temperature

kappa_drop = 2.5 # Kappa drop factor, Similar to drops reported by Baillie et al. (2015) in their opacity calculations for H2O snow-line

#%%
'''
Section to determine the correct R array for use with FARGO, this means writing in 3 ghost zone cells either side of the R limits
and matching the spacing between Radii
'''
#R_FARGO_edges = np.loadtxt('domain_y.dat')
#
#R_FARGO_midpoints = []
#
#for i in range(0, len(R_FARGO_edges) - 1):
#    
#    midpoint = 0.5*(R_FARGO_edges[i] + R_FARGO_edges[i+1])
#    
#    R_FARGO_midpoints.append(midpoint)
#    
#
#R = np.array(R_FARGO_midpoints) #SWITCH (THIS FOR FARGO)


#%%
'''
This section defines the Radius array and functions needed later
'''

R = np.linspace(3, 150, 5e3) #Radius array

alpha_frag_array = np.array([5e-3, 0.01, 0.1]) #This sets the alpha which will be plotted as contours on the final parameter study plot

T_snowline_start = 20 #Inner edge temperature of snowline

dT = 1 #Width of snowline in temperature

Mdot_array = np.logspace(-6, -8, 3)
#Mdot_array = np.array([1.5e-7, 2e-7]) #Array of mass accretion rates to search over [Solar masses per year]


# Opacity (Kappa) function

def kappa_inside(T): # Sets opacity interior to the snowline 
    
    k = (kappa0/kappa_drop)*(T**beta)
    
    return k

def kappa_outside(T): # Sets opacity interior to the snowline 
    
    k = kappa0*(T**beta)
    
    return k


T_snowline_end = T_snowline_start - dT #Outer edge temperature of snowline


def kappa(T): # Sets opacity depending on location in the disc (i.e. whether you are interior, exterior or inside the snowline)


    if T > T_snowline_start: # Interior 
        k = (kappa0/kappa_drop)*(T**beta) 
        
    elif T < T_snowline_end: # Exterior
        k = (kappa0)*(T**beta)
        
        
    else: # Inside snowline
        m = - (kappa_outside(T_snowline_end) - kappa_inside(T_snowline_start))/dT
        c = kappa_outside(T_snowline_end)
        k = m*(T - T_snowline_end) + c
        
    if T > 100: # Sets a maximum constant opacity above a maximum temperature
        k = (kappa0/kappa_drop)*(100**beta)
        
    return k


#%%
# Defining functions


def rec(x): 
    
    recipricol = 1/x
    return recipricol


def Sigma(T, R): #Surface density 
    
    sig = np.sqrt((kb*M)/(mu*G))*(1/(np.pi*Q0))*((R*au)**(-3/2))*np.sqrt(T)
    
    return sig



def f_Mdot(T, R, Mdot): # Function which we set equal to zero and solve numerically later (Essentially radiative balance)
    
    F = ((3/(8*np.pi*sb))*G*M*Mdot*sm*(rec((R*au)**3))) - (((T**4) - (Tirr**4))/((kappa(T)*Sigma(T, R)) + rec((kappa(T)*Sigma(T, R)))))
    
    return F



def T_analytical_Mdot(R, Mdot): # Analytical temperature solution in the optically thick limit with Tirr = 0K, useful for checking low radius temperature

    C = (3/(8*(np.pi**2)))*Mdot*sm*kappa0*(M**(3/2))*np.sqrt((G*kb)/mu)*(1/(sb*Q0))

    Temperature = (C*((R*au)**(-9/2)))**(2/3) 
    
    return Temperature



def T_single_Mdot(R, Mdot): #Returns a value of Temperature for a single Radius (Numerically)
    
    #first_guess = T_analytical_Mdot(R, Mdot)
    first_guess = 50
    
    T_current = sp.optimize.fsolve(f_Mdot, first_guess, args = (R, Mdot), maxfev = int(1e4))
    
    return T_current[0]


def Fj(R, Mdot): #Angular momentum flux
    
    AM_flux = Mdot*sm*np.sqrt(G*M*R*au)
    
    return AM_flux

def tau(T, R): #Optical depth 
    
    t = kappa(T)*Sigma(T, R)
    
    return t

def alpha_gt(T, R): #Gravoturbulent alpha from Rafikov 2015
    Omega = np.sqrt(G*M/((R*au)**3))

    alpha = (8/9)*sb*((np.pi*G*Q0)**6)*((mu/kb)**4)*rec(Sigma(T, R)*kappa(T) + rec(Sigma(T, R)*kappa(T)))*(Sigma(T, R)**5)*rec(Omega**7)*(1 - ((Tirr/T)**4))
    
    return alpha


def alpha_gt_data(Temperature_data, Sigma_data, Tau_data, R): #Gravoturbulent alpha which takes arrays as inputs
    
    Omega = np.sqrt(G*M/((np.array(R)*au)**3))
    
    alpha = (8/9)*sb*((np.pi*G*Q0)**6)*((mu/kb)**4)*rec(np.array(Tau_data) + rec(np.array(Tau_data)))*(np.array(Sigma_data)**5)*rec(Omega**7)*(1 - ((Tirr/np.array(Temperature_data))**4))
    
    return alpha



'''
The following temperature solver returns temperature data for array of R (Radius) for a set Mdot  (Single mass accretion rate)
It Uses fsolve for temperatures above the snowline (+5K and above the snowline) and brent q below that
'''

    
def T_Mdot_combination(R, Mdot): 
    
    T = [[], [], []] # The 3 temperature solutions
    
    T_lowest = Tirr - 1 #For now, becuase it cannot be lower than 12. Maybe later add an fsolve temperature lower
    # especially for water ice line which is at ~150K nowhere near 12K to improve speed. 
    
    
    T_highest = T_snowline_start + 5 # 5K hotter than warm side of snowline
    
    T_guesses = np.linspace(T_lowest, T_highest, 500)
    
    #T_guess = T_analytical_Mdot(R[0], Mdot) #Can cause problems
    T_guess = 500
    
    T_current = T_single_Mdot(R[0], Mdot) + 1 #Sets first T_current so the following if statement works for first R
    
    for R_current in R: #Loop over radius
        
        if T_current > T_highest: #Solve for temperature for +5K above snowline and ABOVE, here there is only 1 solution
            
            T_current = sp.optimize.fsolve(f_Mdot, T_guess, args = (R_current, Mdot), maxfev = int(1e5))[0] # Solve y for single value of x (x_current here)
            # [0] at end of fsolve causes it to return the value
            T_guess = T_current
            
            T[0].append(T_current) # These solutions are the same 
            T[1].append(T_current)
            T[2].append(T_current)
            
        else: #Solve for temperature for +5K above snowline and BELOW, here there is only 1 solution
            
            cross = []
        
            for i in range(0, len(T_guesses) - 1):
                
                
                f_current = f_Mdot(T_guesses[i], R_current, Mdot)
                f_next = f_Mdot(T_guesses[i+1], R_current, Mdot)
                
                mult = f_current*f_next
                
                
                if mult < 0:
            
                    cross.extend((i, i+1))
                    
                    
            for i in range(0, len(cross), 2):
                
                T_solution = sp.optimize.brentq(f_Mdot, T_guesses[cross[i]], T_guesses[cross[i+1]], args = (R_current, Mdot))
                
                if len(cross) == 2: # If the 3 temperature solutions are the same
                    
                        T[0].append(T_solution)
                        T[1].append(T_solution)
                        T[2].append(T_solution)
                        
                elif len(cross) == 4: # If we have 2 identical temperature solutions and 1 different one (We never observed this)
                    
                    if i == 0:
                    
                        T[0].append(T_solution)
                        T[2].append(T_solution)
    
                        
                    else:
                        
                        T[1].append(T_solution)
                        
                else: # For 3 distinct solutions (i.e. in the snowline for viscous instability)
                    
                    if i == 0:
                        
                        T[0].append(T_solution)
                        
                    elif i == 2:
                        
                        T[1].append(T_solution)
                        
                    else:
                        
                        T[2].append(T_solution)
                        
        
    return T  


    

def tau_data(R): #Returns optical depth data 
    
    Tau_data = [[], [], []]
    
    for i in range(len(Tau_data)):
    
        Tau_data_current = []
        
        for j in range(len(R)):
            
            Tau_current = tau(Temperature_data[i][j], R[j])
        
            Tau_data_current.append(Tau_current)
            
        Tau_data[i].append(Tau_data_current)
        
    return Tau_data
