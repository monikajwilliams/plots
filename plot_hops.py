#!/usr/bin/env python
import matplotlib
import os,sys
import numpy as np
import matplotlib.pyplot as plt
import plotHOPS as pH
from numpy import linalg as LA


def msd_hops(nd=1,
	 tstep=50,
	 axis=2,
	 nwaters=100,
	 pDefect=True,
	 intervalps = 5,
	 lenMSDps = 100, 
	 filename2="dOs",
	 filename3="msd_hops",
	 dropps = 10.0,
	 filename="hops"):


	np.set_printoptions(threshold=10000)
	
	infoname = "info.txt"
	fs2ps = 0.001 

	timestep = tstep #Timestep of simulation times timestep of saving to the trajectory
	dropTraj = dropps/(timestep*fs2ps) #in timestep/fram units
	lenMSD = int(lenMSDps/(timestep*fs2ps))
#        interval = int((intervalps/(fs2ps*timestep)))
        #print "Interval = %d steps " % (interval)
	
	# => Trajectory Relevant Data <= #
	hops = np.load(filename)
	#dOs = np.load(filename2)
	print "Loaded Data"

        #total_steps = len(dOs) 
        num_transfers = len(hops)
        interval = len(hops)/4
        num_origins = num_transfers/interval
        
        all_msds = np.zeros(num_transfers)
        normalization = np.zeros(num_transfers)
        count = 0
        for origin in range(0,num_transfers/2,interval):
            msd_seg = []
            dhop = 0.0
            nsamples1 = 0
            nsamples2 = 0
            for m in range(origin,num_transfers-origin):
                count += 1
                nsamples = num_transfers - origin
                if count == 0:
                    nsamples1 += nsamples
                if count == 1:
                    nsamples2 += nsamples
                dhop += hops[origin+m]
                print "=================="
                print hops[origin+m]
                print origin + m
                msd_seg.append(dhop)

            all_msds[:len(msd_seg)] += msd_seg
            normalization[:len(msd_seg)] += np.ones(len(msd_seg))

        #nSDs = len(hops)
	#dSDs = interval
        #print "CHECK"
        #print count            
        #print len(hops)
	##for i in range(count):
	#for i in range(len(hops)):
	#    normalization[:(nSDs-(dSDs*(i)))] += 1
        #    #print nSDs - (dSDs*i)
        print "=========="
        print normalization
        print "=========="
	MSD = np.divide(all_msds,normalization)

        MSD = np.array(MSD)
	MSD.dump(filename3)
	info = "Dropped Trajectory = %d\n " %(dropps)
	info1 = "Time interval between origins = %d\n " %(intervalps)
	hs  = open("info.txt","a")
	hs.write(info)
	hs.write(info1)
	hs.close() 

filenameO = "nhop.pdf"
msd_hops()
D = pH.dcplot(nd=1,
	      tstep=50,
	      fitrange=[0,20],
              nwaters=100,
	     # plotAxis=[0.0,100.0,0.0,100.0],
	      filename2="msd_hops",
              filename=filenameO)
print D



