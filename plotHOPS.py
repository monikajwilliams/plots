#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

def dcplot(nd = 1,
	   fn_in = "msd_vals",
	   fn_out = "msd_hops.png",
	   timestep = 50,
	   fitrange = [5.0,20.0],
	   axis_range = 0,
	   ):	


        # => Loading MSD <= #

	MSD = np.load(fn_in)

	# => Unit Conversions <= #

	nm22a2 = 10.0*10.0
	fs2ps = 0.001

	timesL = np.arange(0.0,len(MSD))
	times = np.array(timesL*timestep*fs2ps)
	
	MSD_A = np.array(MSD*nm22a2)

        # => Fitting Line Segment <= #

	nps = fitrange[0]
        nskip = int(nps/(fs2ps*timestep))
	eps = fitrange[1]
	eskip = int(eps/(fs2ps*timestep))
	
	plt.clf()
	plt.plot(times,MSD_A, alpha=0.3,c='b')

	if eskip-nskip >= len(times):
		m, b = np.polyfit(times, MSD_A, 1)
		print "Skip range is longer than length of simulation!"
		D_value = False
	else:
		m, b = np.polyfit(times[nskip:eskip], MSD_A[nskip:eskip], 1)
		D_value = True

	plt.plot(times[nskip:eskip],(times[nskip:eskip]*m)+b, alpha=1.0,c='r')
	plt.xlabel(r'$\mathrm{t\ (ps)}$',fontsize=20)
	plt.ylabel('$\mathrm{MSD\ O^*\ (\AA^2)}$',fontsize=20)

	MSD_max = max(MSD_A)
	if axis_range != 0:
		plt.axis((axis_range[0],axis_range[1],axis_range[2],axis_range[3]))

	plt.savefig(fn_out) 

	if D_value == False:
		D = 0
		print "Error in simulation, diffusion coefficient is not valid!"
	else:
		D = (m/(2*nd))

	return D
	
def multiplot(MSDs,
	   filename,
           labels,
	   nd=1,
	   tstep=50,
	   fitrange=[20.0,100.0],
	   axis_range=0,
	   ):	
        Ds = [] 
        k = 0
        for MSD in MSDs:

        	#Conversion factors
        	nm22a2 = 10.0*10.0
        	fs2ps = 0.001
        
        	np.set_printoptions(threshold=10000)
        	timestep = tstep
        	
        	
        
        	#used for water!
        	#if nd == 1:
        	#    timesL = np.arange((len(MSD[0,:])))
        	#if nd == 2:
        	#    timesL = np.arange((len(MSD[0,:])))
        	#if nd == 3:
        	#    timesL = np.arange((len(MSD[0,:])))
        
        	timesL = np.arange(len(MSD[0]))
        	times = np.array(timesL*timestep*fs2ps)
        	
        	MSD_A = np.array(MSD*nm22a2)
        
        	#avg_MSD = np.sum(MSD_A, axis=0)/nwaters
        	
        	val1 = 0.0
        	ind1 = 0
        
        	#<========(Line Segment to Fit)========>
        	nps = fitrange[0]
        	nskip = int(nps/(fs2ps*timestep))
        	eps = fitrange[1]
        	eskip = int(eps/(fs2ps*timestep))
        	
        	#plt.plot(times,MSD_A[0], alpha=0.3,c='b', label = labels[k])
        	plt.plot(times,MSD_A[0],label = labels[k])
        
        
        	if eskip-nskip >= len(times):
        		m, b = np.polyfit(times, MSD_A[0], 1)
        		print "Skip range is longer than length of simulation!"
        		D_value = False
        	else:
        		m, b = np.polyfit(times[nskip:eskip], MSD_A[0,nskip:eskip], 1)
        		D_value = True
        
        
        	#plt.plot(times[nskip:eskip],(times[nskip:eskip]*m)+b, alpha=1.0,c='r')
        
        	plt.xlabel(r'$\mathrm{t\ (ps)}$', fontsize=20)
        	plt.ylabel(r'$\mathrm{MSD\ O^*\ (\AA^2)}$',fontsize=20)
        	MSD_max = max(MSD_A[0])
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
        #plt.axis([0.0,100.0,0.0,1000.0])
	plt.savefig(filename) 
        return Ds

def ldMSD(filename):
	MSD = np.load(filename)
        
        return MSD

def transfers(fn_transfers = "transfers",
              fn_plot = "transfers.png",
              timestep = 50.0,
             ):

    # => Parameters <= #
    avg_range = 100.0 # ps
    axis = 2
   
    # => Unit Conversions <= #
 
    fs2ps = 0.001
    ps2steps = 1.0/(timestep*fs2ps)    
    step2ps = 1.0/ps2steps

    transfers = np.load(fn_transfers)

    traj_len = len(transfers)*timestep*fs2ps
    transfers_ps =  np.sum(transfers)/traj_len

    return transfers_ps
  #  times = np.linspace(0.0,traj_len,len(transfers))

  #  # => Calculating Moving Average of Transfer Rate <= #

  #  p = int(avg_range*ps2steps/(2.0))
  #  y = np.zeros(len(transfers)) # to be transfer rate
  #  transfers = transfers.astype(float)

  #  for n,trans in enumerate(transfers):
  #      if n < p:
  #          y[n] = trans/(avg_range)
  #      if n >= len(transfers)-p:
  #          y[n] = trans/(avg_range)
  #      else:
  #          y[n] = sum(transfers[n-p:n+p])/(avg_range)

  #  Z = np.load("Ostar")[:,axis]
  #  Z -= min(Z)

  #  # This step accounts for differences in drop_traj
  #  positions  = Z[np.shape(Z)[0]-np.shape(times)[0]:]

  #  # => Plotting position against transfer rate <= #

  #  plt.clf()
  #  fig = plt.figure()
  #  plt.subplot(211)
  #  plt.plot(times[p:-p],positions[p:-p])
  #  plt.ylabel(r"$\mathrm{Position\ (\AA)}$",fontsize=16)
  # 
  #  plt.subplot(212) 
  #  plt.plot(times[p:-p],y[p:-p],'-')
  #  plt.xlabel(r"$\mathrm{Time\ (ps)}$",fontsize=16)
  #  plt.ylabel(r"$\mathrm{Transfer\ Rate\ (/ps)}$",fontsize=16)
  #  plt.savefig(fn_plot)
  #  
  #  # => Calculating Correlation Between Velocity and Transfer Rate <= #
  #      
  #  n = p
  #  #smooth_positions = np.zeros(len(positions)-2.0*n)
  #  smooth_positions = np.zeros(len(positions))
  #  times_vel = np.linspace(0.0,traj_len,len(smooth_positions))
  #  times_pos = np.linspace(0.0,traj_len,len(positions))

  #  print "Smoothing position data"
  #  for i,val in enumerate(positions):
  #      if i > n and i < len(positions)-2.0*n:
  #          smooth_positions[i] = np.mean(positions[i-n:i+n])
  #      else:
  #          smooth_positions[i] = val

  #  dy = step2ps
  #  print "Taking gradient of position data"
  #  velocities = np.gradient(smooth_positions[:],dy)

  #  p = 2*p
  #  #inds_V = np.argsort(y[p:-p])
  #  velocities2 = abs(velocities[p:-p])
  #  inds_V = np.argsort(velocities2)
  #  times_vel2 = times_vel[p:-p]
  #  
  #  plt.clf()
  #  plt.figure()
  #  plt.plot(velocities2[inds_V[:-1]],y[inds_V[:-1]],'o',alpha=0.2)
  #  plt.xlabel(r"$\mathrm{Velocity\ (\AA/ps)}$")
  #  plt.ylabel(r"$\mathrm{Transfer\ Rate\ (/ps)}$")
# #   plt.axis([0.0,0.2,5.0,12.0])
  #  #plt.plot(y2[inds_V],velocities2[inds_V],'-')
  #  #plt.axis([5.0,10.0,-0.2,0.2])
  #  plt.savefig("Aug26_corr_test.png")

  #  plt.clf()
  #  plt.subplot(211)
  #  plt.plot(times_pos,smooth_positions)
  #  plt.plot(times,positions,c='b',alpha=0.2)

  #  plt.subplot(212)
  #  plt.plot(times[2*p:-2*p],velocities[2*p:-2*p],'-')
  #  plt.savefig("vel_test.png")

    return
    

