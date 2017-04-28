#!/usr/bin/env python
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import newcolors as nc
import numpy as np


def contour_plot(filename,
	   filename2="tau",	
	   filename3="angles",	
           nbins = 50,
           ):


	np.set_printoptions(threshold=10000)

        #Loading and symmetrizing data	
	tau2 = np.load(filename2)
        tau3 = tau2[:]*-1.0
	tau = np.append(tau2[:],tau3)
	angle2 = np.load(filename3)
	angle = np.append(angle2[:],angle2[:])

        xedges = np.arange(-0.3,0.3,0.01)
        yedges = np.arange(110.0,200.0,5.0)
        nbins = [xedges,yedges]
        histVals,xedges,yedges = np.histogram2d(tau,angle,nbins)

        xc = 0.5 * (xedges[:-1] + xedges[1:]) 
        yc = 0.5 * (yedges[:-1] + yedges[1:]) 

        E = histVals
        X, Y = np.meshgrid(xc,yc, indexing='ij')
     
        fig = plt.figure() 
        interval = 0.1*np.max(E)
        levels = np.arange(0.0,2500.0,150.0)
        CS = plt.contourf(X,Y,E,cmap=nc.viridis,levels = levels)
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'$\mathrm{Count}$',fontsize=20)

	plt.xlabel(r'$\mathrm{\tau}$',fontsize=25)
	plt.ylabel(r'$\mathrm{O^\dagger-H-O^*\ angle\ (degrees)}$',fontsize=16)
        plt.axis(fontsize=12)
        xmin = -0.25
        xmax =  0.25
        ymin =  120.0
        ymax =  190.0
        plt.axis([xmin,xmax,ymin,ymax])
	plt.savefig(filename) 

	return 
	
