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
ny = 30

#%% INTEGRATION
tEnd = 10.
nb_steps = 5000

theta_flow = 1.;
accuracy_u = 0.0001;
accuracy_v = 0.0001;
accuracy_p = 0.0001;
clean_pressure = False;

theta_th = 1.
accuracy_th = 0.0001


#%% FLUID PROPERTIES
lamb = 0
rho = 1000
cp = 4000
mu = 50


#%% BODY
R = 0.1*ny
body = [[False for j in range(ny)] for i in range(nx)]
for i in range(nx):
	for j in range(ny):
		if (i-0.5*nx)**2 + (j-0.5*ny)**2 < R**2:
			body[i][j] = True


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