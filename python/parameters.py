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
ny = 20

#%% INTEGRATION
tEnd = 10.
nb_steps = 10000

theta_th = 0.5
accuracy_th = 0.0001


#%% FLUID PROPERTIES
lamb = 0
rho = 1000
cp = 4000
mu = 10


#%% BOUNDARY CONDITIONS
T_ext = 0.
#UP-DOWN-LEFT-RIGHT
U_BC = [0. for i in range(nx)] + [0. for i in range(nx)] + [1. for i in range(ny)] + [1. for i in range(ny)]
V_BC = [0. for i in range(nx)] + [0. for i in range(nx)] + [0. for i in range(ny)] + [0. for i in range(ny)]
P_BC = [0. for i in range(nx)] + [0. for i in range(nx)] + [0. for i in range(ny)] + [0. for i in range(ny)]


#%% INITIAL STATE
U0 = np.zeros((nx,ny))
V0 = np.zeros((nx,ny))
P0 = np.zeros((nx,ny))
T0 = np.zeros((nx,ny))

for i in range(nx):
    for j in range(ny):
        #r2 = (i-nx/2.)**2 + (j-ny/2.)**2
        T0[i,j] = np.exp(-0.1*(i-0.2*nx)**2)
        #T0[i,j] = 0.

for i in range(nx):
	for j in range(ny):
		U0[i,j] = 1.
		V0[i,j] = 0.


#%% PLOT
plot_temperature = True
plot_pressure = True
plot_velocity = True

max_frames = 100
max_arrows_x = 20
max_arrows_y = 20