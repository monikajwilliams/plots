#!/usr/bin/env python
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import newcolors as nc
import numpy as np


def contour_plot(
           fn_plot = "EZ_plot.pdf",
	   fn_angle = "EZ_angle",	
	   fn_delta = "max_delta",	
           nbins = 50,
           ):


	np.set_printoptions(threshold=10000)

        # => Loading and symmetrizing data <= #

	delta1 = np.load(fn_delta)
        delta2 = delta1[:]*-1.0
	delta = np.append(delta1[:],delta2)

	angle1 = np.load(fn_angle)
	angle = np.append(angle1[:],angle1[:])

        #xedges = np.arange(-0.3,0.3,0.01)
        #yedges = np.arange(110.0,200.0,5.0)
        #nbins = [xedges,yedges]

        histVals,xedges,yedges = np.histogram2d(delta,angle,nbins)

        xc = 0.5 * (xedges[:-1] + xedges[1:]) 
        yc = 0.5 * (yedges[:-1] + yedges[1:]) 

        E = histVals
        X, Y = np.meshgrid(xc,yc, indexing='ij')
     
        fig = plt.figure() 
        interval = 0.1*np.max(E)

#        levels = np.arange(0.0,2500.0,150.0)
        #CS = plt.contourf(X,Y,E,cmap=nc.viridis,levels = levels)

        CS = plt.contourf(X,Y,E,cmap=nc.viridis)
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'$\mathrm{Count}$',fontsize=20)

	plt.xlabel(r'$\mathrm{\delta}$',fontsize=25)
	plt.ylabel(r'$\mathrm{O^\1-O^*-O^2\ angle\ (degrees)}$',fontsize=16)
        plt.axis(fontsize=12)

        #xmin = -0.25
        #xmax =  0.25
        #ymin =  120.0
        #ymax =  190.0
        #plt.axis([xmin,xmax,ymin,ymax])

	plt.savefig(fn_plot) 

	return 

def test():
    
    	contour_plot()

test()
