#!/usr/bin/env python
import os, sys, re, math
import all_msd as am
import plotMSD as pM
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import trackOstar as tO


filename = "Ostar"
filename2 = "Zs_4D_100W.png"
timestep = 50.0
fs2ps = 0.001 
ps2ns = 0.001
dropps = 10.0 #Time to drop from beginning of simulation (ps)
dropTraj = dropps/(timestep*fs2ps) #in timestep/fram units
indTraj = int(dropTraj) #index for dropTraj

# => Trajectory Relevant Data <= #
Ostars = np.load(filename)
print "Loaded Data"
OstarZs = Ostars
OstarZs -= min(OstarZs[indTraj:,2]) 

# This is in timestep times frequency frames are saved to trajectory
times = np.arange(len(OstarZs)-int(dropTraj))
timesL = np.arange(len(OstarZs[indTraj:,2]))
times = np.array(timesL*timestep*fs2ps)

plt.plot(times,OstarZs[indTraj:,2])
plt.xlabel('t [ps]')
plt.ylabel('Z position O$^*$ [$\AA$]')
plt.axis([0.0,5000.0,0.0,30.0])
plt.tight_layout()
plt.savefig(filename2) 
