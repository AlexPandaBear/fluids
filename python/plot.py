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



#%% PARAMETERS

exec(open("parameters.py").read())

file_name = "../data/test"
max_frames = 100



#%% PLOT

S = vf.SM()
S.loadData(file_name)

T = np.zeros((nx,ny,nb_steps+1))
P = np.zeros((nx,ny,nb_steps+1))
U = np.zeros((nx,ny,nb_steps+1))
V = np.zeros((nx,ny,nb_steps+1))

for t in range(nb_steps+1):
    for i in range(nx):
        for j in range(ny):
            T[i,j,t] = S.getTemperatureAt(t, i, j)
            P[i,j,t] = S.getPressureAt(t, i, j)
            U[i,j,t] = S.getXVelocityAt(t, i, j)
            V[i,j,t] = S.getYVelocityAt(t, i, j)


x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
x, y = np.meshgrid(x, y)

"""
plt.figure()

for i in range(nbx):
    for j in range(nby):
        res = s.computeVelocityAt(x[i], y[j], step, periodicity, width)
        U[j,i] = res[0]
        V[j,i] = res[1]
        M[j,i] = res[2]
        A[j,i] = res[3]

h = min(width/(nbx-1), height/(nby-1))

if showVelocityVectors:
    for i in range(nbx):
        for j in range(nby):
            u_star = U[j,i]/U0
            v_star = V[j,i]/U0
            plt.annotate('', xy=( x[i] + h*u_star , y[j] + h*v_star ), xytext=(x[i],y[j]), arrowprops={'arrowstyle': '->', 'lw': 2, 'color': 'white'}, va='center')

plt.contourf(X, Y, M, alpha=.75, cmap='jet')
c = plt.contour(X, Y, M, colors='black')
plt.clabel(c)

plt.show()
"""

def update_plot(frame_number, zarray, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(x, y, zarray[:,:,frame_number], cmap="magma")


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


plot = [ax.plot_surface(x, y, T[:,:,0], color='0.75', rstride=1, cstride=1)]
ax.set_zlim(0,T0_max)


frames = range(nb_steps+1)
if (nb_steps+1 > max_frames):
    frames = [int(i*nb_steps/max_frames) for i in range(max_frames+1)]

animate = animation.FuncAnimation(fig, update_plot, frames, fargs=(T, plot))
plt.show()

#writer = animation.PillowWriter(fps=15)  
#animate.save("test.gif", writer=writer)