import numpy as np
import _fluids as fl


LBM = fl.LBM()
LBM.buildSpace(0., 1., 5, 0., 1., 5)
LBM.buildTime(0., 1., 5)


rho0 = fl.Field([5, 5])
u0 = fl.Field([5, 5, 2])
e0 = fl.Field([5, 5])

for i in range(5):
	for j in range(5):
		rho0[[i,j]] = 1. + np.exp(-((i-2)**2 + (j-2)**2))
		u0[[i,j,0]] = 1.
		u0[[i,j,1]] = 0.
		e0[[i,j]] = 100.

gamma = 1.4

LBM.setInitialState(rho0, u0, e0, gamma)


nu =1.

col = fl.BGK(nu, gamma)
stream = fl.Stream2D([])
out = fl.Field([])

LBM.simulate(col, stream, out)

""" Test buffer protocol
R = np.asarray(rho0)

for i in range(5):
	print([rho0[[i,j]] for j in range(5)])

print("--------")

for i in range(5):
	print([R[i,j] for j in range(5)])
"""