#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 19:05:56 2020

@author: alexandre
"""

#%% IMPORTS


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import _vf as vf



#%% DATA


Lx = 10
Ly = 10
nx = 100
ny = 100

tEnd = 10
nb_steps = 1000

lamb = 0
rho = 1000
cv = 4000
mu = 0

T0 = np.zeros((nx,ny))
for i in range(nx):
    for j in range(ny):
        r2 = (i-nx/2.)**2 + (j-ny/2.)**2
        T0[i,j] = 150*(np.exp(-0.01*r2) + np.sin(0.1*i)**2)

T0_max = 500

U_max = 1
V_max = 1

U0 = np.zeros((nx,ny))
V0 = np.zeros((nx,ny))
for i in range(nx):
    for j in range(ny):
        U0[i,j] = U_max * np.cos(np.pi*((i/nx)-0.5)) * np.sin(np.pi*((2*j/ny)-1))
        V0[i,j] = -V_max * np.cos(np.pi*((j/ny)-0.5)) * np.sin(np.pi*((2*i/nx)-1))

P0_max = 10**5

P0 = np.zeros((nx,ny))
for i in range(nx):
    for j in range(ny):
        P0[i,j] = P0_max


#%% SIMULATION


S = vf.SM()
S.defineTimeParameters(tEnd, nb_steps)
S.defineGridParameters(Lx, Ly, nx, ny)
S.defineFluidProperties(lamb, rho, cv, mu)
S.defineInitialState(U0, V0, P0, T0)
S.defineBoundaryConditions("flux", 0)
S.launchSimulation()
T = S.getResults()

TEMP = np.zeros((nx,ny,nb_steps+1))
for t in range(nb_steps+1):
    for i in range(nx):
        for j in range(ny):
            TEMP[i,j,t] = T[t][i][j]



#%% PLOT


def update_plot(frame_number, zarray, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(x, y, zarray[:,:,frame_number], cmap="magma")


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
x, y = np.meshgrid(x, y)


plot = [ax.plot_surface(x, y, TEMP[:,:,0], color='0.75', rstride=1, cstride=1)]
ax.set_zlim(0,T0_max)


animate = animation.FuncAnimation(fig, update_plot, nb_steps+1, fargs=(TEMP, plot))
plt.show()

#writer = animation.PillowWriter(fps=15)  
#animate.save("test.gif", writer=writer)