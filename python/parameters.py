#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 19:05:56 2020

@author: alexandre
"""

import numpy as np

#Meshing
Lx = 10
Ly = 10
nx = 50
ny = 50

#Time
tEnd = 5
nb_steps = 500

#Fluid
lamb = 0
rho = 1000
cv = 4000
mu = 0

#Initial state
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
        U0[i,j] = U_max * np.cos(np.pi*((i/nx)-0.5)) * np.sin(np.pi*((2*j/ny)-1))
        V0[i,j] = -V_max * np.cos(np.pi*((j/ny)-0.5)) * np.sin(np.pi*((2*i/nx)-1))

P0_max = 10**5

for i in range(nx):
    for j in range(ny):
        P0[i,j] = P0_max

#Error tolerance
pressure_err = 10.