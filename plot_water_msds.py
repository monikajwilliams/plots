#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from numpy import linalg as LA
import mdtraj as mdt


def dcplot(filename,
	   nd=3,
	   fn_in="msd_vals",
	   timestep=50.0,
	   fitrange=[1.0,5.0],
	   axis_range=0,
           nm2a = False,
           ):

   # => Loading data <= #
   
   msd = np.load(fn_in)
   times_steps = np.arange((len(msd[0,:])))
   
   # = > Unit Conversions <= #
   
   # Conversion factors
   if nm2a == True:
       nm22a2 = 100.0
   else:
       nm22a2 = 1.0
   fs2ps = 0.001
   
   times = np.array(times_steps*timestep*fs2ps)
   msd_a = np.array(msd*nm22a2)
   nmols = np.shape(msd)[0]
   avg_msd = np.sum(msd_a, axis=0)/nmols
   
   # setting up fit range
   nps = fitrange[0]
   nskip = int(nps/(fs2ps*timestep))
   eps = fitrange[1]
   eskip = int(eps/(fs2ps*timestep))
   
   # => Plotting <= #	
   
   plt.clf()
   plt.plot(times,avg_msd, alpha=1.0,c='b')
   for n in range(0,nmols,10):
      plt.plot(times,msd_a[n,:], alpha=0.3,c='g')
   
   if eskip-nskip > len(times):
   	print "Skip range is longer than length of simulation!"
   	m, b = np.polyfit(times, avg_msd[0], 1)
   	D_value = False
   else:
   	m, b = np.polyfit(times[nskip:eskip], avg_msd[nskip:eskip], 1)
   	D_value = True
   
   plt.plot(times[nskip:eskip],(times[nskip:eskip]*m)+b, alpha=1.0,c='r')
   plt.xlabel(r'$\mathrm{t\ (ps)}$',fontsize=20)
   plt.ylabel(r'$\mathrm{MSD\ O^*\ (\AA^2)}$',fontsize=20)
   
   msd_max = max(msd_a[0])
   if axis_range != 0:
   	plt.axis((axis_range[0],axis_range[1],axis_range[2],axis_range[3]))
   plt.savefig(filename) 
   
   if D_value == False:
   	D = 0
   	print "Error in simulation, diffusion coefficient is not valid!"
   else:
   	D = (m/(2*nd))
   
   return D
