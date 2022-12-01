#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:58:44 2022

@author: kev
"""

import numpy as np
import matplotlib.pyplot as plt


def nearest_indx(arr, val):
    indx = (np.abs(arr - val)).argmin()
    return indx

# Parameters

gamma = 5/3
R = 8.31e7
kb = 1.38e-16
mp = 1.67e-24
mu = 2.35*mp
NA = 6.02e23
au = 1.495978707e13
M = 2e33
G = 6.67e-8
sm = (2e33)/(365*24*60**2) #Solar masses per year in g per s
Tirr = 12
kappa0 = 1.41e-5
sb = 5.6704e-5
#Tirr = np.load('Tirr.npy')
#kappa0 = np.load('kappa0.npy')
#sb = np.load('sb.npy')
yr = 365*24*60**2


# Code unit definitions------------------------

# Code units
G_SF = 1.0
MSTAR_SF = 1.0
R0_SF = 1.0
R_MU_SF = 1.0
MU0_SF = 1.0

#CGS units
G_CGS = 6.674e-8
MSTAR_CGS = 1.9891e33
R0_CGS = 1.49597871e13 # 1au right now, perhaps change to radius of our snowline later
R_MU_CGS = 3.463103e7
MU0_CGS = 12.5663706143591


'''
Returns scale free units from input cgs units
'''

def T_SF(T_CGS):
    
    Tsf = T_CGS*(G_SF*MSTAR_SF/(R0_SF*R_MU_SF))/(G_CGS*MSTAR_CGS/(R0_CGS*R_MU_CGS))
    
    return Tsf



def T_CGS(T_SF):
    
    Tcgs = T_SF/((G_SF*MSTAR_SF/(R0_SF*R_MU_SF))/(G_CGS*MSTAR_CGS/(R0_CGS*R_MU_CGS)))
    
    return Tcgs

# Returns scale free k0 from input cgs k0

def k0_SF(k0_CGS):
    
    C = (G_SF*MSTAR_SF/(R0_SF*R_MU_SF))/(G_CGS*MSTAR_CGS/(R0_CGS*R_MU_CGS))
    D = MSTAR_CGS/(R0_CGS**2)
    
    return k0_CGS*D/(C**2)



# Energy density
# Returns scale free energy density from input cgs energy density


def e_SF(e_CGS):
    
    top = (G_SF*(MSTAR_SF**2))/(R0_SF**3)
    
    bottom = (G_CGS*(MSTAR_CGS**2))/(R0_CGS**3)    
    
    return e_CGS*(top/bottom)

# Density 
    
def d_SF(d_CGS):
    
    #1/bla because scale free numerator M/R^2 is just 1
    
    denominator = MSTAR_CGS/(R0_CGS**2)
    
    return d_CGS*(1/denominator) #Illistrative

# Radius
    
def r_SF(r_CGS):
    
    return r_CGS/R0_CGS


    
# NB: ENERGY = ENERGY DENSITY from now on 

def energy(T, Sig):
    
    e = (kb*T*Sig)/((gamma - 1)*mu)
    
    return e


def pressure(energy, gamma):
    
    return (gamma - 1)*energy

#Velocity
def v_SF(v_CGS):
    
    v = v_CGS/np.sqrt((G_CGS*MSTAR_CGS)/R0_CGS)
    
    return v
    
'''
Returns cgs units from input scale free temperature
'''  

def d_CGS(d_SF):
    
    #1/bla because scale free numerator M/R^2 is just 1
    
    numerator = MSTAR_CGS/(R0_CGS**2)
    
    return d_SF*numerator #Illistrative


def r_CGS(r_SF):
    
    return r_SF*R0_CGS


def pot_CGS(pot_SF):
    
    return pot_SF*((G_CGS*MSTAR_CGS)/R0_CGS)


def e_CGS(e_SF):
    
    top = (G_CGS*(MSTAR_CGS**2))/(R0_CGS**3)   
    
    bottom = (G_SF*(MSTAR_SF**2))/(R0_SF**3)
    
    return e_SF*(top/bottom)
    

def v_CGS(v_SF):
    
    v = v_SF*np.sqrt((G_CGS*MSTAR_CGS)/R0_CGS)
    
    return v



def H(energy, density, radius, vx):
    
    h = (radius/vx)*np.sqrt(energy*(gamma - 1)/density)
    
    return h


def AR(T, R):
    
    Omega = np.sqrt(G*M/((R*au)**3))
    sound_speed = np.sqrt(kb*T/mu)
    
    return sound_speed/(R*Omega*au)

def rec(x):
    
    return 1/x


def alpha_gt(T, R):
    Omega = np.sqrt(G*M/((R*au)**3))

    alpha = (8/9)*sb*((np.pi*G*Q0)**6)*((mu/kb)**4)*rec(Sigma(T, R)*kappa(T) + rec(Sigma(T, R)*kappa(T)))*(Sigma(T, R)**5)*rec(Omega**7)*(1 - ((Tirr/T)**4))
    
    return alpha

