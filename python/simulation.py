#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 19:05:56 2020

@author: alexandre
"""

import numpy as np
import _vf as vf

exec(open("parameters.py").read())

V = np.sqrt(np.amax(U0**2 + V0**2))
Re = rho * V * max(Lx,Ly) / mu
print("Re = {}".format(Re))

S = vf.SM()
S.defineTimeParameters(tEnd, nb_steps)
S.defineGridParameters(Lx, Ly, nx, ny)
S.defineFluidProperties(lamb, rho, cp, mu)
S.defineBody(body)
S.defineInitialState(U0, V0, T0)
S.defineDynamicBoundaryConditions(U_BC, V_BC, P_BC)
S.defineThermalBoundaryConditions("temp", T_ext)
S.defineFlowIntegrationParameters(theta_flow, accuracy_u, accuracy_v, accuracy_p, clean_pressure)
S.defineThermalIntegrationParameters(theta_th, accuracy_th)
S.launchSimulation()
S.saveData(file_name)