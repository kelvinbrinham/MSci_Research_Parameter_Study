#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 11:29:27 2022

@author: Kelvin Brinham
"""

'''
This parameter study searches for viscous instability over a range of mass accretion rates and radii for steady circumstellar discs.

It consists of 3 scripts, this is part c

Part c finds viscous instability and plots various graphs

Part c can be used (after reloading script a) without script b because it reloads data from files generated by b. In other words one 
can manipulate data produced previously in this script.

The first few sections just plot temperature etc. if you want the main viscous instability search purpose of this script run it from 
--
Fj vs Sigma plots 

Finding viscous instability
--
Onwards.


NB: This parameter study uses cgs units unless stated otherwise
NB: Interior/Exterior means at a lower/higher radius than the snowline
'''


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize

import Parameter_study_post_a

from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes



# Loading data
G = np.load('G.npy')
M = np.load('M.npy')
Mdot_array = np.load('Mdot_array.npy')
Q0 = np.load('Q0.npy')
R = np.load('R.npy')
Radii_inner = np.load('Radii_inner.npy')
Radii_inner_index = np.load('Radii_inner_index.npy')
Radii_middle = np.load('Radii_middle.npy')
Radii_outer = np.load('Radii_outer.npy')
Radii_outer_index = np.load('Radii_outer_index.npy')
Radii_snowline = np.load('Radii_snowline.npy')
Radius_array = np.load('Radius_array.npy')
Sigma_data = np.load('Sigma_data.npy')
T = np.load('T.npy')
T_snowline_end = np.load('T_snowline_end.npy')
Tau_data = np.load('Tau_data.npy')
Temperature_data = np.load('Temperature_data.npy')
Temperature_data_inner = np.load('Temperature_data_inner.npy')
Temperature_data_outer = np.load('Temperature_data_outer.npy')
alpha_data = np.load('alpha_data.npy')
alpha_frag = np.load('alpha_frag.npy')
au = np.load('au.npy')
beta = np.load('beta.npy')
dT = np.load('dT.npy')
kappa0 = np.load('kappa0.npy')
kappa_drop = np.load('kappa_drop.npy')
kb = np.load('kb.npy')
mu = np.load('mu.npy')
sb = np.load('sb.npy')
sm = np.load('sm.npy')
#gamma = np.load('gamma.npy')
Tirr = np.load('Tirr.npy')

# Colour blind colours for graphs
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']



#%%
'''
Calculate position of snowline for a continuum of Mdot (GRAPHICALLY)
This section is more for visualisation

Plots Temperature vs R for continuum of Mdots
Plots Sigma vs R for continuum of Mdots
    
'''
#Plot Temperature

for j in range(3):
    
#    plt.figure() #Remove to plot all solutions on 1 plot
    for i in range(0, len(Mdot_array), 100): #Plot every (last number in range) Mdot #PROPER
#    for i in range(0, len(Mdot_array), len(Mdot_array) +1): #Plot every (last number in range) Mdot #REPORT PLOTS

        
        plt.loglog(R, T[i][j], label = '$\dot{M}$ = ' + str(Mdot_array[i]) + '$M_{\odot} yr^{-1}$')
        plt.xlabel('R /au')
        plt.ylabel('T')
        plt.xlim(3, 150)
        plt.ylim(10, 500)
        plt.title('Temperature vs Radius, solution ' + str(j+1))
#        plt.loglog(Radii_inner[i], 100, '|', markersize = 1e4)
#        plt.loglog(Radii_outer[i], 100, '|', markersize = 1e4)
        plt.legend()
#%% Making report temperature plot

fig,ax = plt.subplots()#figsize=(10,10))
ax.set_xlabel('Radius [au]', fontsize = 12)
ax.set_ylabel('Temperature [K]', fontsize = 12)     
ax.annotate('Q$_0$ = ' + str(Q0) + ', T$_{irr}$ = ' + str(Tirr) + 'K' + '\n' + '$\dot{M}$ = ' + str(Mdot_array[0]) + '$M_{\odot} yr^{-1}$', fontsize = 12, xy=(0.03, 0.85), xycoords='axes fraction', bbox = dict(facecolor = 'white', edgecolor = 'black', alpha = 1))
#ax.grid()
    
#axins1 = zoomed_inset_axes(ax, zoom = 6, loc = 1)
MARKER = ['x', '+', 'o']

for j in range(0, 3):

    ax.loglog(R, T[0][j], label = 'Solution ' + str(j + 1))
#    axins1.loglog(R, T[0][j])
    ax.legend(loc = 3, markerscale = 4, framealpha = 1)

# SPECIFY THE LIMITS
x1, x2, y1, y2 = 22, 33, 17, 22
axins1.set_xlim(x1, x2) 
axins1.set_ylim(y1, y2)

axins1.axes.xaxis.set_visible(False)
axins1.axes.yaxis.set_visible(False)

#mark_inset(ax, axins1, loc1=3, loc2=4, fc='none', ec="0")
#plt.savefig('T_vs_R_REPORT.png', dpi = 800)


#%%
#------------------------------------------------------------------  
#Plot Surface density
#plt.figure()
plt.xlabel('R')
plt.ylabel('$\Sigma$')
plt.title('Surface Density vs Radius')
plt.annotate('Q$_0$ = ' + str(Q0) + ', T$_{irr}$ = ' + str(Tirr) + 'K', fontsize = 12, xy=(0.05, 0.1), xycoords='axes fraction', bbox = dict(facecolor = 'none', edgecolor = 'black', alpha = 0.5))

for i in range(0, len(Mdot_array), 600):
    
    for j in range(3):
        plt.loglog(R, Sigma_data[i][j], label = '$\dot{M}$ = ' + str(Mdot_array[i]) + '$M_{\odot} yr^{-1}$' + 'Solution ' + str(j + 1))
  
#    plt.figure('Sigma solution 1')
#    plt.loglog(R, Sigma_data[i][1], label = '$\dot{M}$ = ' + str(Mdot_array[i]) + '$M_{\odot} yr^{-1}$' + 'Solution 1')
#    plt.annotate('Q$_0$ = ' + str(Q0) + ', T$_{irr}$ = ' + str(Tirr) + 'K', fontsize = 12, xy=(0.05, 0.1), xycoords='axes fraction', bbox = dict(facecolor = 'none', edgecolor = 'black', alpha = 0.5))
#    plt.legend()
#    
    
#%%
#---------------------------------------------------------------------------
#Plot Optical depth
#plt.figure()
plt.annotate('Q$_0$ = ' + str(Q0) + ', T$_{irr}$ = ' + str(Tirr) + 'K', fontsize = 12, xy=(0.05, 0.05), xycoords='axes fraction', bbox = dict(facecolor = 'none', edgecolor = 'black', alpha = 0.5))
plt.xlabel('R')
plt.ylabel('Optical Depth')
plt.title('Optical Depth vs Radius')
plt.figure('Tau solution 1')
plt.figure('Tau solution 2')
plt.figure('Tau solution 3')


for i in range(len(Mdot_array)):

    # Plot multiple solutions on same graph
    #for j in range(3):
     #   plt.loglog(R, Tau_data[i][j][0], label = '$\dot{M}$ = ' + str(Mdot_array[i]) + '$M_{\odot} yr^{-1}$' + 'Solution ' + str(j + 1))
    
    # Plot mulptiple solutions on seperate graphs
    #plt.figure('Tau solution 1')
    #plt.loglog(R, Tau_data[i][0][0], label = '$\dot{M}$ = ' + str(Mdot_array[i]) + '$M_{\odot} yr^{-1}$' + 'Solution 1')
    
    plt.figure('Tau solution 2')
    plt.loglog(R, Tau_data[i][1][0], label = '$\dot{M}$ = ' + str(Mdot_array[i]) + '$M_{\odot} yr^{-1}$' + 'Solution 2')
    
    #plt.figure('Tau solution 3')
    #plt.loglog(R, Tau_data[i][2][0], label = '$\dot{M}$ = ' + str(Mdot_array[i]) + '$M_{\odot} yr^{-1}$' + 'Solution 3')

        
    #plt.loglog(Radii_inner[i], Tau_data[i][0][0][np.where(R == Radii_inner[i])[0][0]], '|', markersize = 100)#, label = '$R_{inner}$')
    #plt.loglog(Radii_outer[i], Tau_data[i][0][0][np.where(R == Radii_outer[i])[0][0]], '|', markersize = 100)#, label = '$R_{outer}$')

#plt.figure('Tau solution 1')
#plt.loglog(R, np.repeat(1, len(R)), '--', color = 'black', label = '$\u03C4$ = 1')
#plt.legend()

plt.figure('Tau solution 2')
plt.loglog(R, np.repeat(1, len(R)), '--', color = 'black', label = '$\u03C4$ = 1')
plt.legend()

#plt.figure('Tau solution 3')
#plt.loglog(R, np.repeat(3, len(R)), '--', color = 'black', label = '$\u03C4$ = 1')
#plt.legend()
#%%
#------------------------------------------------------------------
#Plot Alpha
#plt.figure()
plt.annotate('Q$_0$ = ' + str(Q0) + ', T$_{irr}$ = ' + str(Tirr) + 'K', fontsize = 12, xy=(0.05, 0.9), xycoords='axes fraction', bbox = dict(facecolor = 'none', edgecolor = 'black', alpha = 0.5))
plt.xlabel('R')
plt.ylabel(r'$\alpha_{gt}$')
plt.title(r'$\alpha_{gt}$ vs Radius')
plt.figure('alpha solution 1')
plt.figure('alpha solution 2')
plt.figure('alpha solution 3')


for i in range(len(Mdot_array)):
    
    # Plot multiple solutions on same graph
    #for j in range(3):
        
   #     plt.loglog(R, alpha_data[i][j], label = '$\dot{M}$ = ' + str(Mdot_array[i]) + '$M_{\odot} yr^{-1}$' + 'Solution ' + str(j + 1))
        
    # Plot mulptiple solutions on seperate graphs
    plt.figure('alpha solution 1')
    plt.loglog(R, alpha_data[i][0], label = '$\dot{M}$ = ' + str(Mdot_array[i]) + '$M_{\odot} yr^{-1}$' + 'Solution 1')
    
    plt.figure('alpha solution 2')
    plt.loglog(R, alpha_data[i][1], label = '$\dot{M}$ = ' + str(Mdot_array[i]) + '$M_{\odot} yr^{-1}$' + 'Solution 2')
    
    plt.figure('alpha solution 3')
    plt.loglog(R, alpha_data[i][2], label = '$\dot{M}$ = ' + str(Mdot_array[i]) + '$M_{\odot} yr^{-1}$' + 'Solution 3')

        
    #plt.loglog(Radii_inner[i], Sigma_data[i][0][np.where(R == Radii_inner[i])[0][0]], '|', markersize = 100)#, label = '$R_{inner}$')
    #plt.loglog(Radii_outer[i], Sigma_data[i][0][np.where(R == Radii_outer[i])[0][0]], '|', markersize = 100)#, label = '$R_{outer}$')

plt.figure('alpha solution 1')
plt.loglog(R, np.repeat(0.1, len(R)), '--', color = 'black', label = r'$\alpha_{gt, frag}$ = ' + str(alpha_frag))
plt.legend()

plt.figure('alpha solution 2')
plt.loglog(R, np.repeat(0.1, len(R)), '--', color = 'black', label = r'$\alpha_{gt, frag}$ = ' + str(alpha_frag))
plt.legend()

plt.figure('alpha solution 3')
plt.loglog(R, np.repeat(0.1, len(R)), '--', color = 'black', label = r'$\alpha_{gt, frag}$ = ' + str(alpha_frag))
plt.legend()

#%%
'''
Fj vs Sigma plots 

Finding viscous instability

'''


Fj_vs_Sigma_grad = [] #Gradient in Fj vs Sigma plot, of which there is one for every Radius in Radius_array, each with 3 solutions
#[[[], [], []], [[], [], []]...] i.e. [[], [], []] for each Radius in Radius_array for each solution, each point corresponts to an Mdot


#fig,ax = plt.subplots()#figsize=(10,10))
#ax.set_xlabel('Surface Density (Constant Radius) [gcm$^{-2}$]', fontsize = 12)
#ax.set_ylabel('Angular Momentum Flux [gcm$^2$s$^{-2}$]', fontsize = 12)

#ax.set_xticks(ticks=y)
#ax.set_xticklabels(labels=ylabs)




for j in range(0, len(Radius_array)): #Loop over snowline radii
#for j in range(0, 500, 100): 
    
    Fj_data = []
    Sigma_data_single_Radius = [[], [], []]
    
    Fj_vs_Sigma_single_Radius_grad = [] #[[], [], []] Gradient of Fj vs Sigma for each solution, each value in each list corresponds to an Mdot

    
    for i in range(len(Mdot_array)): # loop over Mdot
        
        # First calculate the Fj at each Mdot
        #Fj_current = Fj(Radius, Mdot_array[i])
        Fj_current = Fj(Radius_array[j], Mdot_array[i])

        Fj_data.append(Fj_current) 
        
        # Second Calculate Sigma at each Mdot
        
        sigma_1 = Sigma_data[i][0][Radii_inner_index[-1] + j]
        sigma_2 = Sigma_data[i][1][Radii_inner_index[-1] + j]
        sigma_3 = Sigma_data[i][2][Radii_inner_index[-1] + j]
        
        
        Sigma_data_single_Radius[0].append(sigma_1)
        Sigma_data_single_Radius[1].append(sigma_2)
        Sigma_data_single_Radius[2].append(sigma_3)
        
        
    for k in range(3):
        
        Fj_vs_Sigma_single_soln_single_Radius_grad = np.gradient(Fj_data, Sigma_data_single_Radius[k])
        Fj_vs_Sigma_single_Radius_grad.append(Fj_vs_Sigma_single_soln_single_Radius_grad)
    
    
    Fj_vs_Sigma_grad.append(Fj_vs_Sigma_single_Radius_grad)
        
'''
The following commented code is used to plot Fj vs sigma, it is commented because the code currently just finds viscous instability
which is then represented on the final Mdot vs R plot which is produced at the end of this script.
'''

#    if j == 200:
#            
#
#        Sigma_data_single_Radius_insetdata = Sigma_data_single_Radius[1]
#        Fj_data_insetdata = Fj_data
#
#  
#    ax.loglog(Sigma_data_single_Radius[1], Fj_data, label = 'Radius = ' + '{:.0f}'.format(R[Radii_inner_index[-1] + j]) + 'au')
#    ax.legend(loc = 1, framealpha = 1)
#    
#    
#axins1 = zoomed_inset_axes(ax, zoom = 3.7, loc = 3)
##axins1.loglog(Sigma_data_single_Radius[1], Fj_data, label = 'Radius = ' + '{:.0f}'.format(R[Radii_inner_index[-1] + j]) + 'au', color = 'C2')
#axins1.loglog(Sigma_data_single_Radius_insetdata, np.asarray(Fj_data) - 0.33e39, label = 'Radius = ' + '{:.0f}'.format(R[Radii_inner_index[-1] + j]) + 'au', color = 'C2')
#
#
## SPECIFY THE LIMITS
##x1, x2, y1, y2 = 63, 95, 2.7e39, 6.4e39 
#x1, x2, y1, y2 = 142, 180, 1e39, 1.8e39
#axins1.set_xlim(x1, x2) 
#axins1.set_ylim(y1, y2)
#axins1.axes.xaxis.set_visible(False)
#axins1.axes.yaxis.set_visible(False)
#
#        
## IF SET TO TRUE, TICKS ALONG 
## THE TWO AXIS WILL BE VISIBLE
##plt.xticks(visible=False)
##plt.yticks(visible=False)
#
#mark_inset(ax, axins1, loc1=2, loc2=4, fc='none', ec="0")
#plt.savefig('Fj_vs_Sigma_REPORT.png', dpi = 800)
#
#
#    plt.figure('Fj vs Sigma solution 1')
#    plt.plot(Sigma_data_single_Radius[0], Fj_data, label = 'Radius = ' + str(R[Radii_inner_index[-1] + j]) + 'au')
#    plt.title('$F_j$ vs $\Sigma$ solution 1')
#    plt.xlabel('$\Sigma$')
#    plt.ylabel('$F_j$')
#    plt.legend()
#    
#    plt.figure('Fj vs Sigma solution 2')
#
#
#    plt.ylim(1e38, 1e42)
#    plt.loglog(Sigma_data_single_Radius[1], Fj_data, label = 'Radius = ' + '{:.0f}'.format(R[Radii_inner_index[-1] + j]) + 'au')
#    plt.title('Angular Momentum Flux vs $\Sigma$ solution 2')
#    plt.xlabel('Surface Density (Constant Radius) [gcm$^{-2}$]', fontsize = 15)
#    plt.ylabel('$F_j$')
#    plt.ylabel('Angular Momentum Flux [gcm$^2$s$^{-2}$]', fontsize = 15)
#    plt.legend()
#    plt.annotate('Q$_0$ = ' + str(Q0) + ', T$_{irr}$ = ' + str(Tirr) + 'K, ' + '$\kappa_0$ = 1.41x10$^{-5}$$cm^{2}g^{-1}K^{-2}$', fontsize = 10, xy=(0.05, 0.05), xycoords='axes fraction', bbox = dict(facecolor = 'white', edgecolor = 'black', alpha = 1))
#    plt.legend(loc = 1, framealpha = 1)
#    plt.annotate('Decreasing Radius',fontsize = 15,
#            xy=(0.85, 0.95), xycoords='axes fraction',
#            xytext=(0.65, 0.95), textcoords='axes fraction',
#            arrowprops=dict(facecolor='black', shrink=0.05),
#            horizontalalignment='right', verticalalignment='center')
#
#    plt.savefig('Fj_vs_Sigma_example_for_viva_presentation.png', dpi = 800)
#    
#
#    
#    plt.figure('Fj vs Sigma solution 3')
#    plt.plot(Sigma_data_single_Radius[2], Fj_data, label = 'Radius = ' + str(R[Radii_inner_index[-1] + j]) + 'au')
#    plt.title('$F_j$ vs $\Sigma$ solution 3')
#    plt.xlabel('$\Sigma$')
#    plt.ylabel('$F_j$')
#    plt.legend()

#    
#    for k in range(3):
#        
#        Fj_vs_Sigma_single_soln_single_Radius_grad = np.gradient(Fj_data, Sigma_data_single_Radius[k])
#        Fj_vs_Sigma_single_Radius_grad.append(Fj_vs_Sigma_single_soln_single_Radius_grad)
#    
#    
#    Fj_vs_Sigma_grad.append(Fj_vs_Sigma_single_Radius_grad)
            



'''
Finding negative gradients ==> Instability
'''

Mdot_vs_R_data = [[], [], []] # each [] for a solution, each one is s.t [R, Mdot, #(see below)]

# 0 = Stable, Doesnt Fragment
# 1 = Unstable, Doesnt Fragment
# 2 = Stable, Fragments
# 3 = Unstable, Fragments



for i in range(len(Radius_array)): #Loop over each Radius
    
    for j in range(3): #Loop over each solution
        
        for k in range(len(Mdot_array)): #Loop over each Mdot
        
            
            if Fj_vs_Sigma_grad[i][j][k] > 0 and alpha_data[k][j][Radii_inner_index[-1] + i] < alpha_frag: #Stable, Doesnt Fragment
            
                Mdot_vs_R_data[j].extend([[Radius_array[i], Mdot_array[k], 0]])
                
            
            elif Fj_vs_Sigma_grad[i][j][k] < 0 and alpha_data[k][j][Radii_inner_index[-1] + i] < alpha_frag: #Unstable, Doesnt Fragment
                
                Mdot_vs_R_data[j].extend([[Radius_array[i], Mdot_array[k], 1]])
                
            
            elif Fj_vs_Sigma_grad[i][j][k] > 0 and alpha_data[k][j][Radii_inner_index[-1] + i] > alpha_frag: #Stable, Fragments
                
                Mdot_vs_R_data[j].extend([[Radius_array[i], Mdot_array[k], 2]])
                
            
            else:
                
                Mdot_vs_R_data[j].extend([[Radius_array[i], Mdot_array[k], 3]]) #Unstable, Fragments
                


'''
Tau = 1 line data

Used for plotting optical depth = 1 line on final plot

'''        
Tau_1_data = [[[], []], [[], []], [[], []]] #[Mdot, R] for each Tau = 1 solution at each Mdot, i.e. the 3 tau solutions for solution 2 of temperature


for i in range(len(Mdot_array)):
    
    Tau_modified = np.array(Tau_data[i][1][0]) - 1
    
    R_midpoint_data = []
    
    for j in range(len(R) - 1):
        
        Tau_modified_current = Tau_modified[j]
        Tau_modified_next = Tau_modified[j+1]
        
        mult = Tau_modified_current*Tau_modified_next
        
        if mult < 0:
            
            R_midpoint = 0.5*(R[j+1] + R[j])
            R_midpoint_data.append(R_midpoint)
            
    
        
    if len(R_midpoint_data) == 1:
        
        Tau_1_data[0][1].append(R_midpoint_data[0])
        Tau_1_data[1][1].append(R_midpoint_data[0])
        Tau_1_data[2][1].append(R_midpoint_data[0])
        
        Tau_1_data[0][0].append(Mdot_array[i])
        Tau_1_data[1][0].append(Mdot_array[i])
        Tau_1_data[2][0].append(Mdot_array[i])
        
    elif len(R_midpoint_data) == 2:
        
        Tau_1_data[0][1].append(R_midpoint_data[0])
        Tau_1_data[1][1].append(R_midpoint_data[1])
        Tau_1_data[2][1].append(R_midpoint_data[0])
        
        Tau_1_data[0][0].append(Mdot_array[i])
        Tau_1_data[1][0].append(Mdot_array[i])
        Tau_1_data[2][0].append(Mdot_array[i])
        
    elif len(R_midpoint_data) == 3:

        
        Tau_1_data[0][1].append(R_midpoint_data[0])
        Tau_1_data[1][1].append(R_midpoint_data[1])
        Tau_1_data[2][1].append(R_midpoint_data[2])
        
        Tau_1_data[0][0].append(Mdot_array[i])
        Tau_1_data[1][0].append(Mdot_array[i])
        Tau_1_data[2][0].append(Mdot_array[i])
        
    else:
        
        print('More than 3 solutions!')
        pass
    



'''

alpha = alpha_frag line data 

Used for plotting fragmentation boundary on final plot.

The alpha_frag lines which are plotted can be set in script a in alpha_frag_array 
'''

def alpha_frag_line_data(alpha_fragment):
    

    alpha_alpha_frag_data = [[[], []], [[], []], [[], []]] #[Mdot, R] where alpha = alpha_frag for each solution for each alpha_frag
    
    R_midpoint_data_length = [] #Length of R_midpoint for each Mdot
    
    for i in range(len(Mdot_array)):
        
        alpha_modified = np.array(alpha_data[i][1]) - alpha_fragment #[i] for which Mdot then [1] for 2nd temperature solution
        
        
        R_midpoint_data = []
        
        for j in range(len(R) - 1):
            
            alpha_modified_current = alpha_modified[j]
            alpha_modified_next = alpha_modified[j+1]
            
            
            mult = alpha_modified_current*alpha_modified_next
            
            if mult < 0:
                
                R_midpoint = 0.5*(R[j+1] + R[j])
                R_midpoint_data.append(R_midpoint)
             
                
        R_midpoint_data_length.append(len(R_midpoint_data))
        
            
        if len(R_midpoint_data) == 1:
            
            
            alpha_alpha_frag_data[0][1].append(R_midpoint_data[0])
            alpha_alpha_frag_data[1][1].append(R_midpoint_data[0])
            alpha_alpha_frag_data[2][1].append(R_midpoint_data[0])
            
            alpha_alpha_frag_data[0][0].append(Mdot_array[i])
            alpha_alpha_frag_data[1][0].append(Mdot_array[i])
            alpha_alpha_frag_data[2][0].append(Mdot_array[i])
            
        elif len(R_midpoint_data) == 2:
                    
            
            alpha_alpha_frag_data[0][1].append(R_midpoint_data[0])
            alpha_alpha_frag_data[1][1].append(R_midpoint_data[1])
            alpha_alpha_frag_data[2][1].append(R_midpoint_data[0])
            
            alpha_alpha_frag_data[0][0].append(Mdot_array[i])
            alpha_alpha_frag_data[1][0].append(Mdot_array[i])
            alpha_alpha_frag_data[2][0].append(Mdot_array[i])
            
        elif len(R_midpoint_data) == 3:
    
                    
            alpha_alpha_frag_data[0][1].append(R_midpoint_data[0])
            alpha_alpha_frag_data[1][1].append(R_midpoint_data[1])
            alpha_alpha_frag_data[2][1].append(R_midpoint_data[2])
            
            alpha_alpha_frag_data[0][0].append(Mdot_array[i])
            alpha_alpha_frag_data[1][0].append(Mdot_array[i])
            alpha_alpha_frag_data[2][0].append(Mdot_array[i])
            
        else:
            
            pass
    
    '''
    Breaking the alpha_alpha_frag_data into solution interior to, inside and exterior to the snowline
    I do this by splitting it where there are 3 solutions
    '''
    #Finding where to break alpha_alpha_frag_data
    R_midpoint_data_length_number_of_solution_changes_index = []
    
    for i in range(len(R_midpoint_data_length) - 1):
        
        R_midpoint_length_current = R_midpoint_data_length[i]
        R_midpoint_length_next = R_midpoint_data_length[i+1]
        
        if R_midpoint_length_current == 0 or R_midpoint_length_next == 0:
            
            pass
        
        else:
        
            if R_midpoint_length_current != R_midpoint_length_next:
            
                R_midpoint_data_length_number_of_solution_changes_index.append(i)
        
    
    if len(R_midpoint_data_length_number_of_solution_changes_index) > 2:
        
        print('Number of solutions changes more than twice!')
                
    global alpha_alpha_frag_data_interior
    global alpha_alpha_frag_data_exterior
    #Breaking alpha_alpha_frag_data, i just pick 1st alpha_alpha_frag solution for interior and exterior to the snowline because all should be equivilant there
    
    #alpha_alpha_frag_data interior to the snowline, i.e. at lower R than the snowline
    
    alpha_alpha_frag_data_interior = [alpha_alpha_frag_data[0][0][:(R_midpoint_data_length_number_of_solution_changes_index[0] + 1)], alpha_alpha_frag_data[0][1][:(R_midpoint_data_length_number_of_solution_changes_index[0] + 1)]]
    
    #alpha_alpha_frag_data inside the snowline
    #alpha_alpha_frag_data_inside = 
    
    #alpha_alpha_frag_data exterior to the snowline, i.e. at higher R than the snowline
    
    #------------------- BODGE HERE -------------------------------
    #alpha_alpha_frag_data_exterior = [alpha_alpha_frag_data[0][0][(R_midpoint_data_length_number_of_solution_changes_index[1] + 1):], alpha_alpha_frag_data[0][1][(R_midpoint_data_length_number_of_solution_changes_index[1] + 1):]]
    alpha_alpha_frag_data_exterior = [alpha_alpha_frag_data[0][0][(R_midpoint_data_length_number_of_solution_changes_index[1] + 3):], alpha_alpha_frag_data[0][1][(R_midpoint_data_length_number_of_solution_changes_index[1] + 3):]] #NEW BODGE TO AVOID THE 2 outer 
    # dodgy points that give 1 point interior to the snowline



#%%------------------------------------------------------------------------------------------------------------------------
'''

Mdot vs R graph plotting with filled regions (just solution 2 points plotted for now) 

'''
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
                  
fig, ax = plt.subplots()
ax.loglog(Radii_inner, Mdot_array, CB_color_cycle[5])#, label = 'Inner Radius')
ax.loglog(Radii_outer, Mdot_array, CB_color_cycle[1])#, label = 'Outer Radius')
#-----------Arrows for poster-----------
plt.annotate('Snowline',fontsize = 18,
            xy=(45, 3.7e-7), xycoords='data',
            xytext=(0.7, 0.6), textcoords='axes fraction',
            arrowprops=dict(facecolor='black', shrink=0.05),
            horizontalalignment='left', verticalalignment='top')
plt.annotate('Optically thick', fontsize = 12, xy = (10.5, 8.5e-7), textcoords='data')
plt.annotate('Optically thin', fontsize = 12, xy = (24, 8.5e-7), textcoords='data')
#---------------------------------------
#ax.set(xlabel = 'Radius /au', ylabel = '$\dot{M}$ /$M_{\odot} yr^{-1}$') #, title = '$\dot{M}$ vs Radius')
ax.set_ylabel('Mass accretion rate [$M_{\odot} yr^{-1}$]', fontsize = 12)
ax.set_xlabel('Radius [au]', fontsize = 12)

Mdot_vs_R_data_unstable_no_fragment = []
Mdot_vs_R_data_unstable_fragment = []


#Array of inner/outer most radii that is unstable for points that don't/do fragment for each Mdot all within 
#the snowline
# e.g. Radii_inner_unstable_no_fragment goes as inner most radius for each Mdot that is unstable and doesnt 
#fragment within the snowline

Radii_inner_unstable_no_fragment = []
Radii_inner_unstable_fragment = []
Radii_outer_unstable_no_fragment = []
Radii_outer_unstable_fragment = []

Mdot_unstable_no_fragment = []
Mdot_unstable_fragment = []


for i in range(len(Mdot_vs_R_data[1])): #Just solution 2 rn, change this first index to see other solutions
   
    
    if Mdot_vs_R_data[1][i][2] == 1: #Unstable, Doesn't Fragment
        
        #Add all unstable non fragmenting [R, Mdot] lists to a single list
        Mdot_vs_R_data_unstable_no_fragment.append([Mdot_vs_R_data[1][i][0], Mdot_vs_R_data[1][i][1]])
        
        
    elif Mdot_vs_R_data[1][i][2] == 3: #Unstable, Fragments
        
        #Add all unstable fragmenting [R, Mdot] lists to a single list
        Mdot_vs_R_data_unstable_fragment.append([Mdot_vs_R_data[1][i][0], Mdot_vs_R_data[1][i][1]])
        
        
    
# Sort above lists by accending Mdot
Mdot_vs_R_data_unstable_no_fragment = sorted(Mdot_vs_R_data_unstable_no_fragment, key = lambda x: x[1])
Mdot_vs_R_data_unstable_fragment = sorted(Mdot_vs_R_data_unstable_fragment, key = lambda x: x[1])

#Radii_inner_unstable_no_fragment.append(Mdot_vs_R_data_unstable_no_fragment[0][0])
Mdot_vs_R_data_unstable_no_fragment.insert(0, [0, 0])
Mdot_vs_R_data_unstable_fragment.insert(0, [0, 0])
Mdot_vs_R_data_unstable_no_fragment.append([np.inf, np.inf])
Mdot_vs_R_data_unstable_fragment.append([np.inf, np.inf])


for i in range(len(Mdot_vs_R_data_unstable_no_fragment) - 1):
    
    if Mdot_vs_R_data_unstable_no_fragment[i][1] != Mdot_vs_R_data_unstable_no_fragment[i+1][1]:
        
        Radii_inner_unstable_no_fragment.append(Mdot_vs_R_data_unstable_no_fragment[i+1][0])
        Radii_outer_unstable_no_fragment.append(Mdot_vs_R_data_unstable_no_fragment[i][0])
        
        Mdot_unstable_no_fragment.append(Mdot_vs_R_data_unstable_no_fragment[i][1])

        

for i in range(len(Mdot_vs_R_data_unstable_fragment) - 1):
    
    if Mdot_vs_R_data_unstable_fragment[i][1] != Mdot_vs_R_data_unstable_fragment[i+1][1]:
        
        Radii_inner_unstable_fragment.append(Mdot_vs_R_data_unstable_fragment[i+1][0])
        Radii_outer_unstable_fragment.append(Mdot_vs_R_data_unstable_fragment[i][0])  
        
        Mdot_unstable_fragment.append(Mdot_vs_R_data_unstable_fragment[i][1])


        
#Trim the bodged infinity and 0 off lists 
del Radii_inner_unstable_no_fragment[-1]
del Radii_outer_unstable_no_fragment[0]
del Radii_inner_unstable_fragment[-1]
del Radii_outer_unstable_fragment[0]   

del Mdot_unstable_no_fragment[0]
del Mdot_unstable_fragment[0]
        
#Adding first highest Mdot solution of non fragmenting points to fragmenting points to fill in the gap between them
#when filling
#-------------------------------------------------------------------------------

# Finding Mdot 1 lower than the lowest Mdot in Mdot_unstable_fragment
# Index of lowest Mdot in Mdot_unstable_fragment in Mdot_array
Mdot_unstable_fragment_1_lower = Mdot_array[np.where(Mdot_array == Mdot_unstable_fragment[0])[0][0] + 1]

Radii_inner_unstable_fragment.insert(0, Radii_inner_unstable_no_fragment[np.where(Mdot_unstable_no_fragment == Mdot_unstable_fragment_1_lower)[0][0]])
Radii_outer_unstable_fragment.insert(0, Radii_inner_unstable_no_fragment[np.where(Mdot_unstable_no_fragment == Mdot_unstable_fragment_1_lower)[0][0]])

Mdot_unstable_fragment.insert(0, Mdot_unstable_fragment_1_lower)


Mdot_unstable_no_fragment_1_higher = Mdot_array[np.where(Mdot_array == Mdot_unstable_no_fragment[-1])[0][0] - 1]

Radii_inner_unstable_no_fragment.append(Radii_outer_unstable_fragment[np.where(Mdot_unstable_fragment == Mdot_unstable_no_fragment_1_higher)[0][0]])
Radii_outer_unstable_no_fragment.append(Radii_outer_unstable_fragment[np.where(Mdot_unstable_fragment == Mdot_unstable_no_fragment_1_higher)[0][0]])

Mdot_unstable_no_fragment.append(Mdot_unstable_no_fragment_1_higher)


#Filling in final end gap in top right (CAUTION THIS WILL NOT WORK PROPERLY WHEN THERE IS NO FRAGMENTATION)
Radii_outer_unstable_fragment[-1] = Radii_outer[0]

#-------------------------------------------------------------------------------


#Fillingw
ax.fill_betweenx(Mdot_unstable_no_fragment, Radii_inner_unstable_no_fragment, Radii_outer_unstable_no_fragment, alpha = 0.5, color = CB_color_cycle[0], label = 'Unstable and Doesn$\'$t Fragment')
ax.fill_betweenx(Mdot_unstable_fragment, Radii_inner_unstable_fragment, Radii_outer_unstable_fragment, alpha = 0.5, color = CB_color_cycle[7], label = 'Unstable and Fragments')

# Plotting the Tau = 1 lines
ax.loglog(Tau_1_data[0][1], Tau_1_data[0][0], '.', markersize = 0.2, color = 'black')#, label = 'Optical depth = 1')
ax.loglog(Tau_1_data[1][1], Tau_1_data[1][0], '.', markersize = 0.2, color = 'black')
ax.loglog(Tau_1_data[2][1], Tau_1_data[2][0], '.', markersize = 0.2, color = 'black')

#Full alpha = alpha_frag solution 
#ax.loglog(alpha_alpha_frag_data[0][1], alpha_alpha_frag_data[0][0], '--', color = 'black', label = r'$\alpha = \alpha_{frag}$ = ' + str(alpha_frag))
#ax.loglog(alpha_alpha_frag_data[1][1], alpha_alpha_frag_data[1][0], '--', color = 'black')
#ax.loglog(alpha_alpha_frag_data[2][1], alpha_alpha_frag_data[2][0], '--', color = 'black')

#alpha = alpha_frag solution just interior and exterior to the snowline 


#alpha_frag_array = np.array([5e-3, 5e-2, 0.01, 0.1])
alpha_frag_array = np.array([0.1])

colors = ['black', 'green', 'grey', 'orange']

for k in range(len(alpha_frag_array)):
    
    alpha_frag = alpha_frag_array[k]
    alpha_frag_line_data(alpha_frag)

    
    ax.loglog(alpha_alpha_frag_data_interior[1], alpha_alpha_frag_data_interior[0], '--', color = colors[k], label = 'Fragmentation boundary')#r'$\alpha_{frag}$ = ' + str("{:.2f}".format(alpha_frag)))
    ax.loglog(alpha_alpha_frag_data_exterior[1], alpha_alpha_frag_data_exterior[0], '--', color = colors[k]) 


ax.annotate('Q$_0$ = ' + str(Q0) + ', T$_{irr}$ = ' + str(Tirr) + 'K, ' + '$\kappa_0$ = 1.41x10$^{-5}$$cm^{2}g^{-1}K^{-2}$', fontsize = 10, xy=(0.3, 0.3), xycoords='axes fraction', bbox = dict(facecolor = 'none', edgecolor = 'black', alpha = 1))
ax.legend(loc = 'lower right', prop={'size': 12}) #10 is default font size



#fig.savefig('Mdot_vs_R_REPORT.png', dpi = 800)

#%%
'''
Save arrays from this parameter study as files to be opened in Converting_code_units_1.0 and to be used later as initial conditions
in FARGO

'''

# Density, middle solution here (i.e. the solution exhibiting viscous instability)


Density = np.asarray(Sigma_data[0][1]) # 1st [] gives the Mdot (Although only 1 here so pick [0]), 2nd [] gives the solution
np.save('Density', Density)

Temperature = np.asarray(T[0][1])
np.save('Temperature', Temperature)

Radius = np.asarray(R)
np.save('Radius', Radius)


'''
Getting t_cool as a function of Radius to compare with final data (This is plotted there)
'''

t_cool = (kb/mu)*Density*(Tau_data[0][1][0] + 1/Tau_data[0][1][0])/(2*(gamma - 1)*sb*((Temperature**2) + (Tirr**2))*(Temperature + Tirr))

np.save('t_cool', t_cool)




