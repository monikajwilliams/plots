#!/usr/bin/env python
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import statsmodels.api as sm
from numpy import linalg as LA
import mdtraj as mdt

np.set_printoptions(threshold=10000)
r = np.arange(10.0)

a2i = -0.2281 #kcal/(mol*A^2)
a4i = 1.09    #kcal/(mol*A^4)
a6i = 0.2341  #kcal/(mol*A^6)
a8i = 0.3254  #kcal/(mol*A^8)

a2 = a2i*cal2J*(10**2) #kJ/(mol*A^2)
a4 = a4i*cal2J*(10**4) #kJ/(mol*A^4)
a6 = a6i*cal2J*(10**6) #kJ/(mol*A^6)
a8 = a8i*cal2J*(10**8) #kJ/(mol*A^8)


V = a2*(r**2))+a4*(r**4)+a6*(r**6))+a8*(r**8)
plt.plot(r,V)
plt.xlabel(r'Radius ($\AA$)')
plt.ylabel('Potential (kJ/mol)')
plt.savefig('potential.pdf')
