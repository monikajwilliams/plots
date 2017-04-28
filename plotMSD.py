#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from numpy import linalg as LA
import mdtraj as mdt

def dcplot(filename,
	   nd=1,
	   fn_in="msd_vals",
	   timestep=50.0,
	   fitrange=[1.5,4.0],
	   axis_range=0,
           nm2a = True,
           ):

	
        # => Loading data <= #

        msd = np.load(fn_in)
	times_steps = np.arange(len(msd[0]))

        # = > Unit Conversions <= #

	# Conversion factors
        if nm2a == True:
	    nm22a2 = 100.0
        else:
	    nm22a2 = 1.0
	fs2ps = 0.001

	times = np.array(times_steps*timestep*fs2ps)
	msd_a = np.array(msd*nm22a2)

        # setting up fit range
	nps = fitrange[0]
	nskip = int(nps/(fs2ps*timestep))
	eps = fitrange[1]
	eskip = int(eps/(fs2ps*timestep))

        # => Plotting <= #	

	plt.clf()
	plt.plot(times,msd_a[0], alpha=0.3,c='b')

	if eskip-nskip > len(times):
		print "Skip range is longer than length of simulation!"
		m, b = np.polyfit(times, msd_a[0], 1)
		D_value = False
	else:
		m, b = np.polyfit(times[nskip:eskip], msd_a[0,nskip:eskip], 1)
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
	
def multiplot(msds,
	   filename,
           labels,
	   nd=1,
	   tstep=50,
	   fitrange=[20.0,100.0],
	   axis_range=0,
	   ):	
        Ds = [] 
        k = 0
        for msd in msds:

        	#Conversion factors
        	nm22a2 = 10.0*10.0
        	fs2ps = 0.001
        
        	np.set_printoptions(threshold=10000)
        	timestep = tstep
        	
        	
        
        	#used for water!
        	#if nd == 1:
        	#    times_steps = np.arange((len(msd[0,:])))
        	#if nd == 2:
        	#    times_steps = np.arange((len(msd[0,:])))
        	#if nd == 3:
        	#    times_steps = np.arange((len(msd[0,:])))
        
        	times_steps = np.arange(len(msd[0]))
        	times = np.array(times_steps*timestep*fs2ps)
        	
        	msd_a = np.array(msd*nm22a2)
        
        	#avg_msd = np.sum(msd_a, axis=0)/nwaters
        	
        	val1 = 0.0
        	ind1 = 0
        
        	#<========(Line Segment to Fit)========>
        	nps = fitrange[0]
        	nskip = int(nps/(fs2ps*timestep))
        	print "nskip = %d" %(nskip)
        	eps = fitrange[1]
        	eskip = int(eps/(fs2ps*timestep))
        	print "eskip = %d" %(eskip)
        	
        	#plt.plot(times,msd_a[0], alpha=0.3,c='b', label = labels[k])
        	plt.plot(times,msd_a[0],label = labels[k])
        
        
        	if eskip-nskip >= len(times):
        		m, b = np.polyfit(times, msd_a[0], 1)
        		print "Skip range is longer than length of simulation!"
        		D_value = False
        	else:
        		m, b = np.polyfit(times[nskip:eskip], msd_a[0,nskip:eskip], 1)
        		D_value = True
        
        
        	#plt.plot(times[nskip:eskip],(times[nskip:eskip]*m)+b, alpha=1.0,c='r')
        
        	plt.xlabel(r'$\mathrm{t\ (ps)}$',fontsize=20)
        	plt.ylabel(r'$\mathrm{msd\ O^*\ (\AA^2)}$',fontsize=20)
        	msd_max = max(msd_a[0])
        	if axis_range != 0:
        		plt.axis((axis_range[0],axis_range[1],axis_range[2],axis_range[3]))
        
        	if D_value == False:
        		D = 0
        		print "Error in simulation, diffusion coefficient is not valid!"
        	else:
        		D = (m/(2*nd))
                        Ds.append(D)
                k += 1
	plt.tight_layout() 
        plt.legend(loc=4)
        plt.axis([0.0,100.0,0.0,1000.0])
	plt.savefig(filename) 
        return Ds

