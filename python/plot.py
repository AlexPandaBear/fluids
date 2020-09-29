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
import _vf as vf

exec(open("parameters.py").read())



#%% DATA LOADING

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

X = np.zeros((nx, ny))
Y = np.zeros((nx, ny))

for i in range(nx):
    for j in range(ny):
        X[i,j] = x[i]
        Y[i,j] = y[j]



#%% VELOCITY PLOT

if plot_velocity:
    M = np.zeros((nx,ny,nb_steps+1))
    arr_x = np.zeros((nx,ny,nb_steps+1))
    arr_y = np.zeros((nx,ny,nb_steps+1))

    h = min(Lx/(nx-1), Ly/(ny-1))
    u0_max = np.amax(U0)+1
    v0_max = np.amax(V0)+1

    coef_u = h/u0_max
    coef_v = h/v0_max

    for i in range(nx):
        for j in range(ny):
            for t in range(nb_steps+1):
                M[i,j,t] = np.sqrt(U[i,j,t]**2 + V[i,j,t]**2)
                arr_x[i,j,t] = X[i,j] + coef_u*U[i,j,t]
                arr_y[i,j,t] = Y[i,j] + coef_v*V[i,j,t]


    frames = range(nb_steps+1)
    if (nb_steps+1 > max_frames):
        frames = [int(i*nb_steps/max_frames) for i in range(max_frames+1)]

    arrows_x = range(nx)
    if (nx > max_arrows_x):
        arrows_x = [int(i*(nx-1)/max_arrows_x) for i in range(max_arrows_x+1)]

    arrows_y = range(ny)
    if (ny > max_arrows_y):
        arrows_y = [int(i*(ny-1)/max_arrows_y) for i in range(max_arrows_y+1)]


    fig = plt.figure()
    ax = fig.add_subplot(111)

    def update_velocity_plot(frame_number):
        ax.clear()
        ax.contourf(X, Y, M[:,:,frame_number], alpha=.75, cmap='jet')
        c = ax.contour(X, Y, M[:,:,frame_number], colors='black')
        plt.clabel(c)
        for i in arrows_x:
            for j in arrows_y:
                plt.annotate('', xy=(arr_x[i,j,frame_number],arr_y[i,j,frame_number]), xytext=(X[i,j],Y[i,j]), arrowprops={'arrowstyle': '->', 'lw': 2, 'color': 'white'}, va='center')

    animate = animation.FuncAnimation(fig, update_velocity_plot, frames)
    plt.show()




def update_plot(frame_number, zarray, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(X, Y, zarray[:,:,frame_number], cmap="magma")



#%% PRESSURE PLOT

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