#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import statsmodels.api as sm
from numpy import linalg as LA

density_data_ = {
     1  :  3.2,
     2  :  3.4,
     3  :  3.5,
     4  :  3.6,
     5  :  3.7,
     6  :  3.8,
     7  :  3.9,
     8  :  4.0,
     9  :  4.1,
    10  :  4.2,
    11  :  4.3,
    12  :  4.4,
    13  :  4.5,
    } 

data = np.loadtxt('final_summary.txt',delimiter='|',skiprows=1)

#<====(Plots by Number of Water Molecules)=====>
plt.clf()
#Cut 1&6th densities! Fix these!
for n in range(4,10):
    nD = np.where(data == n+1)
    density = (1.0/ density_data_[n+1])*10
    if len(nD) > 0:
        plt.plot(data[nD[0],0],data[nD[0],2],'o-', label= r'%1.1f $\AA$' %(density))
plt.xlabel('Number of Water Molecules')
plt.ylabel('D O$^*$ [$\AA^2/ps$]')
plt.axis((30,100,0,20))
plt.legend(loc=2)
plt.tight_layout() #Cleans figure makes tight and neat :)
plt.savefig('waters_D.pdf') 

#<====(Plots by Density)=====>
plt.clf()
for n in range(30,110,10):
    nD = np.where(data == n)
    densities1 = data[nD[0],1]
    densities = []
    for density1 in densities1:
        #density = density_data_[density1]
        density = (1.0/ density_data_[density1])*10
        densities.append(density)
    np.array(densities)
    if len(densities) > 0:
         plt.plot(densities,data[nD[0],2],'o-', label='%d waters' %(n))
plt.xlabel(r'Average O-O distance $\AA$')
plt.ylabel(r'D O$^*$ [$\AA^2/ps$]')
plt.legend()
#plt.axis((3.5,4.3,0,20))
plt.tight_layout() #Cleans figure makes tight and neat :)
plt.savefig('Density_D.pdf') 
