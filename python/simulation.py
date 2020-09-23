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
S.defineFluidProperties(lamb, rho, cv, mu)
S.defineInitialState(U0, V0, P0, T0)
S.defineBoundaryConditions("flux", 0)
S.defineErrorTolerance(pressure_err)
S.launchSimulation()
S.saveData("../data/test")