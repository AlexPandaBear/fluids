#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 19:05:56 2020

@author: alexandre
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import _vf as vf

exec(open("parameters.py").read())



#%% DATA LOADING

print("Loading data...")
S = vf.SM()
S.loadData(file_name)

print("Extracting primary fields...")
T = [S.getTemperatureFieldAt(t) for t in range(nb_steps+1)]
P = [S.getPressureFieldAt(t) for t in range(nb_steps+1)]
U = [S.getXVelocityFieldAt(t) for t in range(nb_steps+1)]
V = [S.getYVelocityFieldAt(t) for t in range(nb_steps+1)]

print(type(T[0]))


print("Creating meshgrid...")
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)

X = np.zeros((nx, ny))
Y = np.zeros((nx, ny))

for i in range(nx):
    for j in range(ny):
        X[i,j] = x[i]
        Y[i,j] = y[j]



#%% VELOCITY PLOT

print("Plotting velocity field...")
if plot_velocity:
    if (nb_steps+1 > max_frames):
        frames = [int(i*nb_steps/max_frames) for i in range(max_frames+1)]
    else:
        frames = range(nb_steps+1)

    if (nx > max_arrows_x):
        arrows_x = [int(i*(nx-1)/max_arrows_x) for i in range(max_arrows_x+1)]
    else:
        arrows_x = range(nx)
    
    if (ny > max_arrows_y):
        arrows_y = [int(i*(ny-1)/max_arrows_y) for i in range(max_arrows_y+1)]
    else:
        arrows_y = range(ny)

    M = np.zeros((nx,ny,nb_steps+1))
    arr_x = np.zeros((nx,ny,nb_steps+1))
    arr_y = np.zeros((nx,ny,nb_steps+1))

    h = min(Lx/(max_arrows_x-1), Ly/(max_arrows_y-1))
    u0_max = np.sqrt(np.amax(U0**2 + V0**2))

    coef_u = h/u0_max
    coef_v = h/u0_max

    for i in range(nx):
        for j in range(ny):
            for t in range(nb_steps+1):
                M[i,j,t] = np.sqrt(U[i,j,t]**2 + V[i,j,t]**2)
                arr_x[i,j,t] = X[i,j] + coef_u*U[i,j,t]
                arr_y[i,j,t] = Y[i,j] + coef_v*V[i,j,t]


    fig = plt.figure()
    ax = fig.add_subplot(111)
    div = make_axes_locatable(ax)
    cax = div.append_axes('right', '5%', '5%')

    def update_velocity_plot(frame_number):
        ax.clear()
        cax.cla()
        
        cf = ax.contourf(X, Y, M[:,:,frame_number], alpha=.75, cmap='jet', vmin=0., vmax=2.)
        c = ax.contour(X, Y, M[:,:,frame_number], colors='black')
        
        plt.colorbar(cf, cax=cax, shrink=0.9)
        plt.clabel(c)
        
        for i in arrows_x:
            for j in arrows_y:
                ax.annotate('', xy=(arr_x[i,j,frame_number],arr_y[i,j,frame_number]), xytext=(X[i,j],Y[i,j]), arrowprops={'arrowstyle': '->', 'lw': 2, 'color': 'white'}, va='center')

    animate = animation.FuncAnimation(fig, update_velocity_plot, frames)
    plt.show()




def update_plot(frame_number, zarray, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(X, Y, zarray[:,:,frame_number], cmap="magma")



#%% PRESSURE PLOT

print("Plotting pressure field...")
if plot_pressure:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


    plot = [ax.plot_surface(X, Y, P[:,:,0], color='0.75', rstride=1, cstride=1)]
    #ax.set_zlim(0,T0_max)


    frames = range(nb_steps+1)
    if (nb_steps+1 > max_frames):
        frames = [int(i*nb_steps/max_frames) for i in range(max_frames+1)]

    animate = animation.FuncAnimation(fig, update_plot, frames, fargs=(P, plot))
    plt.show()



#%% TEMPERATURE PLOT

print("Plotting temperature field...")
if plot_temperature:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    plot = [ax.plot_surface(X, Y, T[:,:,0], color='0.75', rstride=1, cstride=1)]
    ax.set_zlim(0,5)


    frames = range(nb_steps+1)
    if (nb_steps+1 > max_frames):
        frames = [int(i*nb_steps/max_frames) for i in range(max_frames+1)]

    animate = animation.FuncAnimation(fig, update_plot, frames, fargs=(T, plot))
    plt.show()

    #writer = animation.PillowWriter(fps=15)  
    #animate.save("test.gif", writer=writer)