#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 19:05:56 2020

@author: alexandre
"""

import _vf as vf

exec(open("parameters.py").read())

S = vf.SM()
S.defineTimeParameters(tEnd, nb_steps)
S.defineGridParameters(Lx, Ly, nx, ny)
S.defineFluidProperties(lamb, rho, cp, mu)
S.defineInitialState(U0, V0, T0)
S.defineDynamicBoundaryConditions(U_BC, V_BC, P_BC)
S.defineThermalBoundaryConditions("temp", T_ext)
S.defineThermalIntegrationParameters(theta_th, accuracy_th)
S.launchSimulation()
S.saveData(file_name)