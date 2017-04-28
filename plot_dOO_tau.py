#!/usr/bin/env python
import os,sys,math
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import newcolors as nc
import numpy as np


def contour_plot(filename,
	   filename2="tau",	
	   filename3="dOO",	
           plotAxis = 0.0,
           nbins = 50,
           ):

	#Conversion factors
	nm2a = 10.0

	np.set_printoptions(threshold=10000)
	
        #Loading and symmetrizing data	
	tau2 = np.load(filename2)
        tau3 = tau2[:]*-1.0
	tau = np.append(tau2[:],tau3)
	dOO2 = np.load(filename3)
        dOO_nm = np.append(dOO2[:,0,0],dOO2[:,0,0])
        dOO = dOO_nm*nm2a

        histVals,xedges,yedges = np.histogram2d(tau,dOO,nbins,normed=True)
    
        xc = 0.5 * (xedges[:-1] + xedges[1:]) 
        yc = 0.5 * (yedges[:-1] + yedges[1:]) 

        E = histVals
        X, Y = np.meshgrid(xc,yc, indexing='ij')
     
        fig = plt.figure() 
        levels = np.arange(0.0,52.0,2.0)
        CS = plt.contourf(X,Y,E,cmap=nc.viridis,levels=levels)
        #CS = plt.contourf(X,Y,E,cmap=nc.viridis,levels = np.arange(0,np.max(E)+5.0))
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'$\mathrm{Count}$',fontsize=20)

	plt.xlabel(r'$\mathrm{\tau}$',fontsize=25)
	#plt.ylabel(r'$\mathrm{O-O\ distance\ (\AA)}$',fontsize=20)
	plt.ylabel(r'$\mathrm{R_{O^\dagger-O^*}\ (\AA)}$',fontsize=25)
        plt.axis(fontsize=12)
        xmin = -0.25
        xmax = 0.25
        ymin = 2.4
        ymax = 2.8
        plt.axis([xmin,xmax,ymin,ymax])
	plt.savefig(filename) 

        #CS = plt.contourf(X,Y,E,cmap=nc.viridis,levels = np.arange(0,np.max(E)+5.0))
        #levels = np.arange(5.0,50.0,5.0)
        #cbar = plt.colorbar(CS)
        #cbar.ax.set_ylabel('Count')

	#plt.xlabel(r'$\tau$',fontsize=20)
	#plt.ylabel(r'O-O distance [$\AA$]',fontsize=20)
        #plt.axis(fontsize=12)
        #xmin = np.min(X)
        #xmax =  np.max(X)
        #ymin = np.min(Y)
        #ymax = np.max(Y)
        #plt.axis([xmin,xmax,ymin,ymax])
	#plt.savefig(filename) 

	return 
	
def contour_dff_dOO(dir1,
           dir2,
           dir3,
           filename,
	   filename2="tau",	
	   filename3="dOO",	
           plotAxis = 0.0,
           nbins = 50,
           ):

	#Conversion factors
	nm2a = 10.0

	np.set_printoptions(threshold=10000)
    
        os.chdir(dir1)
        #Loading and symmetrizing data	
	d1tau2 = np.load(filename2)
        d1tau3 = d1tau2[:]*-1.0
	d1tau = np.append(d1tau2[:],d1tau3)
	d1dOO2 = np.load(filename3)
        d1dOO_nm = np.append(d1dOO2[:,0,0],d1dOO2[:,0,0])
        d1dOO = d1dOO_nm*nm2a

        xedges1 = np.arange(-0.3,0.3,0.01)
        yedges1 = np.arange(2.30,2.90,0.05)
        nbins1 = [xedges1,yedges1]
        d1histVals,xedges1,yedges1 = np.histogram2d(d1tau,d1dOO,nbins1,normed=True)
        xc1 = 0.5 * (xedges1[:-1] + xedges1[1:]) 
        yc1 = 0.5 * (yedges1[:-1] + yedges1[1:]) 
        X1, Y1 = np.meshgrid(xc1,yc1, indexing='ij')
    
        os.chdir(dir2)
	d2tau2 = np.load(filename2)
        d2tau3 = d2tau2[:]*-1.0
	d2tau = np.append(d2tau2[:],d2tau3)
	d2dOO2 = np.load(filename3)
        d2dOO_nm = np.append(d2dOO2[:,0,0],d2dOO2[:,0,0])
        d2dOO = d2dOO_nm*nm2a

        xedges2 = np.arange(-0.3,0.3,0.01)
        yedges2 = np.arange(2.30,2.90,0.05)
        nbins2 = [xedges1,yedges1]
        d2histVals,xedges2,yedges2 = np.histogram2d(d2tau,d2dOO,nbins2,normed=True)
        xc2 = 0.5 * (xedges2[:-1] + xedges2[1:]) 
        yc2 = 0.5 * (yedges2[:-1] + yedges2[1:]) 
        X2, Y2 = np.meshgrid(xc2,yc2, indexing='ij')
        
        E = d2histVals-d1histVals
        X1, Y1 = np.meshgrid(xc1,yc1, indexing='ij')
        os.chdir(dir3) 
        levels = np.arange(-12.0,12.0,0.5)
        levels2 = np.arange(0.0,52.0,2.0)
        CS = plt.contourf(X2,Y2,E,cmap=nc.inferno,levels=levels)

        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'$\mathrm{Count}$',fontsize=20)

	plt.xlabel(r'$\mathrm{\tau}$',fontsize=25)
	plt.ylabel(r'$\mathrm{R_{O^\dagger-O^*}\ (\AA)}$',fontsize=25)
        plt.axis(fontsize=12)
        xmin = -0.25
        xmax = 0.25
        ymin = 2.4
        ymax = 2.8
        plt.axis([xmin,xmax,ymin,ymax])
	plt.savefig(filename) 

def contour_dff_angle(dir1,
           dir2,
           dir3,
           filename,
	   filename2="tau",	
	   filename3="angles",	
           plotAxis = 0.0,
           nbins = 50,
           ):

	#Conversion factors
	nm2a = 10.0

	np.set_printoptions(threshold=10000)
    
        os.chdir(dir1)
        #Loading and symmetrizing data	
	d1tau2 = np.load(filename2)
        d1tau3 = d1tau2[:]*-1.0
	d1tau = np.append(d1tau2[:],d1tau3)
	d1angle2 = np.load(filename3)
        d1angle = np.append(d1angle2,d1angle2)

        xedges1 = np.arange(-0.3,0.3,0.01)
        yedges1 = np.arange(110.0,200.0,5.0)
        nbins1 = [xedges1,yedges1]
        d1histVals,xedges,yedges = np.histogram2d(d1tau,d1angle,nbins1,normed=True)
    
        os.chdir(dir2)
	d2tau2 = np.load(filename2)
        d2tau3 = d2tau2[:]*-1.0
	d2tau = np.append(d2tau2[:],d2tau3)
	d2angle2 = np.load(filename3)
        d2angle = np.append(d2angle2,d2angle2)

        xedges2 = np.arange(-0.3,0.3,0.01)
        yedges2 = np.arange(110.0,200.0,5.0)
        nbins2 = [xedges2,yedges2]
        d2histVals,xedges,yedges = np.histogram2d(d2tau,d2angle,nbins2,normed=True)
        xc = 0.5 * (xedges[:-1] + xedges[1:]) 
        yc = 0.5 * (yedges[:-1] + yedges[1:]) 
        
        E = d2histVals-d1histVals
        X, Y = np.meshgrid(xc,yc, indexing='ij')
    
        os.chdir(dir3) 
        fig = plt.figure() 
        levels = np.arange(-5.0,5.0,1.0)
        CS = plt.contourf(X,Y,E,cmap=nc.magma,levels=levels)
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'$\mathrm{Count}$',fontsize=20)

	plt.xlabel(r'$\mathrm{\tau}$',fontsize=20)
	plt.ylabel(r'$\mathrm{O^\dagger-H-O^*\ angle\ (degrees)}$',fontsize=16)
        plt.axis(fontsize=12)
        xmin = -0.25
        xmax = 0.25
        ymin = 120.0
        ymax = 190.0
        plt.axis([xmin,xmax,ymin,ymax])
	plt.savefig(filename) 
