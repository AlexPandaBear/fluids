#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 19:05:56 2020

@author: alexandre
"""

import numpy as np


file_name = "../data/test"


#%% MESHING
Lx = 20
Ly = 5
nx = 50
ny = 50


#%% INTEGRATION
tEnd = 1
nb_steps = 500


#%% FLUID PROPERTIES
lamb = 1
rho = 1000
cv = 4000
mu = 100


#%% BOUNDARY CONDITIONS
T_ext = 0.
#UP-DOWN-LEFT-RIGHT
U_BC = [0. for i in range(nx)] + [0. for i in range(nx)] + [.3 for i in range(ny)] + [.3 for i in range(ny)]
V_BC = [0. for i in range(nx)] + [0. for i in range(nx)] + [0. for i in range(ny)] + [0. for i in range(ny)]
P_BC = [0. for i in range(nx)] + [0. for i in range(nx)] + [0. for i in range(ny)] + [0. for i in range(ny)]


#%% INITIAL STATE
U0 = np.zeros((nx,ny))
V0 = np.zeros((nx,ny))
P0 = np.zeros((nx,ny))
T0 = np.zeros((nx,ny))

for i in range(nx):
    for j in range(ny):
        r2 = (i-nx/2.)**2 + (j-ny/2.)**2
        T0[i,j] = 150*(np.exp(-0.05*r2) + np.sin(0.1*i)**2)

for i in range(nx):
	for j in range(ny):
		U0[i,j] = .3
		V0[i,j] = 0.


#%% PLOT
plot_temperature = True
plot_pressure = True
plot_velocity = True

max_frames = 200
max_arrows_x = 20
max_arrows_y = 20