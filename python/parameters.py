#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 19:05:56 2020

@author: alexandre
"""

import numpy as np


file_name = "../data/test"


#%% MESHING
Lx = 10
Ly = 10
nx = 100
ny = 100


#%% INTEGRATION
tEnd = 10
nb_steps = 200
pressure_err = 10.


#%% FLUID PROPERTIES
lamb = 1
rho = 1000
cv = 4000
mu = 1


#%% INITIAL STATE
U0 = np.zeros((nx,ny))
V0 = np.zeros((nx,ny))
P0 = np.zeros((nx,ny))
T0 = np.zeros((nx,ny))

for i in range(nx):
    for j in range(ny):
        r2 = (i-nx/2.)**2 + (j-ny/2.)**2
        T0[i,j] = 150*(np.exp(-0.05*r2) + np.sin(0.1*i)**2)

T0_max = 500

U_max = 1
V_max = 1

for i in range(nx):
	for j in range(ny):
		i_star = 2.*i/(nx-1) - 1.
		j_star = 2.*j/(ny-1) - 1.
		U0[i,j] = - U_max * np.sin(np.pi*j_star) * np.cos(0.5*np.pi*i_star)
		V0[i,j] = V_max * np.sin(np.pi*i_star) * np.cos(0.5*np.pi*j_star)

P0_max = 10**5

for i in range(nx):
    for j in range(ny):
        P0[i,j] = P0_max


#%% PLOT
plot_temperature = False
plot_pressure = False
plot_velocity = True

max_frames = 100
max_arrows_x = 20
max_arrows_y = 20