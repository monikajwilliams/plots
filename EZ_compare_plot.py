#!/usr/bin/env python
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import newcolors as nc
import numpy as np
import os, sys, math


def contour_delta(
           dir_active,
           dir1,
           dir2,
           fn_plot = "EZ_compare_delta.pdf",
	   fn_angle = "EZ_angle",	
	   fn_max_delta = "max_delta",	
	   fn_min_delta = "min_delta",	
           nbins = 150,
           ):


	np.set_printoptions(threshold=10000)

        # => Loading and symmetrizing data <= #

        os.chdir(dir1)
	min_delta1 = np.load(fn_min_delta)
	max_delta1 = np.load(fn_max_delta)
        delta1 = min_delta1/max_delta1
        #delta1 = min_delta1

	angle1 = np.load(fn_angle)

        xedges = np.arange(-0.1,1.1,0.01)
        yedges = np.arange(50.0,140.0,5.0)
        nbins = [xedges,yedges]

        histVals1,xedges,yedges = np.histogram2d(delta1,angle1,nbins)

        os.chdir(dir2)
	min_delta2 = np.load(fn_min_delta)
	max_delta2 = np.load(fn_max_delta)
        delta2 = min_delta2/max_delta2
        #delta2 = min_delta2

	angle2 = np.load(fn_angle)

        histVals2,xedges,yedges = np.histogram2d(delta2,angle2,nbins)

        xc = 0.5 * (xedges[:-1] + xedges[1:]) 
        yc = 0.5 * (yedges[:-1] + yedges[1:]) 

        E = histVals2-histVals1
        X, Y = np.meshgrid(xc,yc, indexing='ij')
    
        os.chdir(dir_active) 
        plt.clf()
        fig = plt.figure() 
        interval = 0.1*np.max(E)
       
        for n,vals in enumerate(E):
            for m, val in enumerate(vals):
                if abs(val) <= 5.0:
                    E[n,m] = 0.0        

        levels = np.arange(-450.0,490.0,20.0)
        CS = plt.contourf(X,Y,E,cmap=nc.inferno,levels = levels)

        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'$\mathrm{Count}$',fontsize=20)

	plt.xlabel(r'$\mathrm{\delta_{min}/\delta_{max}}$',fontsize=18)
	plt.ylabel(r'$\mathrm{O_1-O^*-O_2\ angle\ (degrees)}$',fontsize=16)
        plt.axis(fontsize=12)

        xmin =  0.0
        xmax =  1.0
        ymin =  60.0
        ymax =  125.0
        plt.axis([xmin,xmax,ymin,ymax])

	plt.savefig(fn_plot) 
        
        xedges2 = np.arange(0.0,1.0,0.05)
        nbins = xedges2
        histdelta1,xedges2 = np.histogram(delta1,nbins)
        histdelta2,xedges2 = np.histogram(delta2,nbins)
        histdelta = histdelta2-histdelta1
        xc2 = 0.5 * (xedges2[:-1] + xedges2[1:]) 
         
        plt.clf()
        plt.figure()

        xmin =  0.0
        xmax =  1.0
        ymin =  min(histdelta) - 100.0
        ymax =  abs(min(histdelta)) + 100.0
        plt.axis([xmin,xmax,ymin,ymax])
        plt.plot(xc2,histdelta,'o-')
	plt.xlabel(r'$\mathrm{\delta_{min}/\delta_{max}}$',fontsize=18)
	plt.ylabel(r'$\mathrm{Bin\ Count}$',fontsize=18)
        plt.savefig("delta_lin_test.pdf")


        x = np.arange(0.0,len(min_delta1),1.0)
        plt.clf()
        plt.figure()
        corr_max1 = delta1[np.argsort(min_delta1)]
        corr_min1 = min_delta1[np.argsort(min_delta1)]*10.0
        corr_max2 = delta2[np.argsort(min_delta2)]
        corr_min2 = min_delta2[np.argsort(min_delta2)]*10.0
        plt.plot(corr_min1, corr_max1,'o',alpha=0.3)
        plt.plot(corr_min2, corr_max2,'o',alpha=0.3)
        #plt.plot(corr_min2, corr_max2-corr_max1,'o')
        
    
        angle1_min = angle1[np.argsort(min_delta1)]
        angle2_min = angle2[np.argsort(min_delta2)]
       # plt.plot(corr_min1, angle1_min,'o')
       # plt.plot(corr_min2, angle2_min,'o')
       # plt.plot(corr_min1, angle2_min-angle1_min,'o')
        
        plt.savefig("min_max_deltas_test.png")

	return 

def contour_dOO(
           dir_active,
           dir1,
           dir2,
           fn_plot = "EZ_compare_dOO.pdf",
	   fn_angle = "EZ_angle",	
	   fn_max_dOO = "max_dOO",	
	   fn_min_dOO = "min_dOO",	
           nbins = 150,
           ):


	np.set_printoptions(threshold=10000)
        nm2a = 10.0

        # => Loading and symmetrizing data <= #

        os.chdir(dir1)
	min_dOO1 = np.load(fn_min_dOO)
	max_dOO1 = np.load(fn_max_dOO)
        dOO1 = min_dOO1/max_dOO1

	angle1 = np.load(fn_angle)

        xedges = np.arange(0.0,1.5,0.01)
        yedges = np.arange(50.0,145.0,5.0)
        nbins = [xedges,yedges]

        histVals1,xedges,yedges = np.histogram2d(dOO1,angle1,nbins)

        os.chdir(dir2)
	min_dOO2 = np.load(fn_min_dOO)
	max_dOO2 = np.load(fn_max_dOO)
        dOO2 = min_dOO2/max_dOO2

	angle2 = np.load(fn_angle)

        histVals2,xedges,yedges = np.histogram2d(dOO2,angle2,nbins)

        xc = 0.5 * (xedges[:-1] + xedges[1:]) 
        yc = 0.5 * (yedges[:-1] + yedges[1:]) 

        E = histVals2-histVals1
        X, Y = np.meshgrid(xc,yc, indexing='ij')
    
        os.chdir(dir_active) 
        plt.clf()
        fig = plt.figure() 
        interval = 0.1*np.max(E)
       
        for n,vals in enumerate(E):
            for m, val in enumerate(vals):
                if abs(val) <= 5.0:
                    E[n,m] = 0.0        

        levels = np.arange(-8950.0,8970.0,550.0)
        CS = plt.contourf(X,Y,E,cmap=nc.inferno,levels = levels)

        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'$\mathrm{Count}$',fontsize=20)

	plt.xlabel(r'$\mathrm{O*-O_{min}\ \AA}$',fontsize=18)
	plt.ylabel(r'$\mathrm{O_1-O^*-O_2\ angle\ (degrees)}$',fontsize=16)
        plt.axis(fontsize=12)

        xmin =  1.0
        xmax =  1.1
        ymin =  70.0
        ymax =  135.0
        plt.axis([xmin,xmax,ymin,ymax])
	plt.savefig(fn_plot) 

        
        xedges2 = np.arange(0.0,1.1,0.05)
        nbins2 = 300
        histdOO1,xedges2 = np.histogram(dOO1,bins=nbins2) 
        histdOO2,xedges2 = np.histogram(dOO2,bins=nbins2) 
        xc2 = 0.5 * (xedges2[:-1] + xedges2[1:]) 
        histvals3 = histdOO2-histdOO1

        plt.clf()
        plt.figure()
        plt.plot(xc2,histvals3,'o-')

        xmin =  0.0
        xmax =  1.1
        ymin =  min(histvals3) - 10.0
        ymax =  max(histvals3) + 10.0
        plt.axis([xmin,xmax,ymin,ymax])
        plt.savefig("dOO_lin_test.pdf")

        yedges2 = np.arange(50.0,155.0,5.0)
        nbins3 = yedges2
        yc2 = 0.5 * (yedges2[:-1] + yedges2[1:]) 
        histangle1,xedges = np.histogram(angle1,bins=nbins3) 
        histangle2,xedges = np.histogram(angle2,bins=nbins3) 
        histvals4 = histangle2-histangle1

        plt.clf()
        plt.figure()
        plt.plot(yc2,histvals4,'o-')

        xmin =  50.0
        xmax =  150.0
        ymin =  min(histvals4) - 100.0
        ymax =  max(histvals4) + 100.0
        plt.axis([xmin,xmax,ymin,ymax])
        plt.savefig("angle_lin_test.pdf")

        plt.clf()
        plt.figure()
        
        corr_dOO1 = min_dOO1[np.argsort(min_dOO1)]*10.0
        corr_angle1 = angle1[np.argsort(min_dOO1)]
        #corr_max1 = dOO1[np.argsort(min_dOO1)]*10.0
        corr_max1 = max_dOO1[np.argsort(min_dOO1)]*10.0

        corr_dOO2 = min_dOO2[np.argsort(min_dOO2)]*10.0
        corr_angle2 = angle2[np.argsort(min_dOO2)]
        #corr_max2 = dOO2[np.argsort(min_dOO2)]*10.0
        corr_max2 = max_dOO2[np.argsort(min_dOO2)]*10.0

        #plt.plot(corr_dOO1,corr_angle1,'o')
        #plt.plot(corr_dOO2,corr_angle2,'o')
        plt.plot(corr_dOO1,corr_max1,'o')
        plt.plot(corr_dOO2,corr_max2,'o')
        plt.savefig("dOO_min_max_test.png")

	return 
