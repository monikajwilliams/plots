#!/usr/bin/env python
import os, sys, re, math
import matplotlib.pyplot as plt
import numpy as np
import time
from datetime import date
from traj_analysis import RDF
from traj_analysis import density as DD

def plot_rdf(
            fn_in,
            waters,
            var_1,
            nd,
            bin_interval = 0.01,
            ones_line = True,
            calculate_distances = False,
            ref_dir = '/Users/monikawilliams/Desktop/current_Research/var_dens/',
            color_ind = 0.06,
            axis_range = [0.0],
            dens_norm = True,
            periodic = False,
            fn_traj='positions.xtc',
            cut = 1E23, # ridiculously high number
            ):

    # => Unit Conversions <= #

    nm2A = 10.0
     
    # Dump files for different types of RDF calculations 
    
    fn_Ostar_O = 'Ostar-O_dists'
    fn_Ostar_H = 'Ostar-H_dists'
    fn_Hstar_O = 'Hstar-O_dists'
    fn_O_O = 'O-O_dists'
    fn_O_H = 'O-H_dists'
    fn_H_H = 'H-H_dists'

    # => Output file with record of options/filenames <= #
    
    date1 = date.today()
    timestamp = "%s-%s" % (str(date1.day), str(date1.month))
    hs = open("RDF.out","a")
   
    if nd == 1: 
        label1 = []
        for i,density in enumerate(var_1):
            avg_OO_dist = (1.0/density)*nm2A
            label_a = r"$\mathrm{Avg.\ O-O\ dist.\ =\ %1.2f\ \AA}$" % (avg_OO_dist)
            label1.append(label_a)
        hs.write("Number of Water Molecules | Density | RDF-Type | Date:\n")
    if nd == 3: 
        label1 = []
        for press in var_1:
            label_a = r"$\mathrm{%1.2f\ GPa}$" % (press)
            label1.append(label_a)
        hs.write("Number of Water Molecules | Pressure | RDF-Type | Date:\n")
    
        
    
    fn_pdf = "RDF_%s_%s.pdf" % (fn_in,timestamp)
    
    fig = plt.figure() 
    
    for n in waters:
        for ind,m in enumerate(var_1):

                if nd == 1: 
        	    os.chdir('%s%dW/%1.1fD/PD/' % (ref_dir,n,m))
        	    os.system("pwd")
                    hs.write(" %d | %1.1f | %s | %s\n" % (n,m,fn_in,timestamp))
                if nd == 3:
                    os.chdir('%spress_%1.2f_GPa/dens/' % (ref_dir,m))
        	    os.system("pwd")
                    hs.write(" %d | %1.1f | %s | %s\n" % (n,m,fn_in,timestamp))

                if calculate_distances == True:

                    if fn_in == fn_Hstar_O:
                        RDF.hstar_o(nwaters = n,
                                    periodic=periodic,
                                    fn_traj=fn_traj,
                                     )
                    elif fn_in == fn_Ostar_O:
                        RDF.ostar_o(nwaters = n,
                                    periodic=periodic,
                                    fn_traj=fn_traj,
                                    )
                    elif fn_in == fn_Ostar_H:
                        RDF.ostar_h(nwaters = n,
                                    periodic=periodic,
                                    fn_traj=fn_traj,
                                    )
                    elif fn_in == fn_O_O:
                        RDF.O_O(nwaters = n,
                                periodic=periodic,
                                fn_traj=fn_traj,
                                cut = cut
                                 )
                    elif fn_in == fn_O_H:
                        RDF.O_H(nwaters = n,
                                periodic=periodic,
                                fn_traj=fn_traj,
                                cut = cut
                                 )
                    elif fn_in == fn_H_H:
                        RDF.H_H(nwaters = n,
                                periodic=periodic,
                                cut = cut,
                                fn_traj=fn_traj,
                                 )

                # Loading Data and Converting from nm to Angstrom

                print "Loading Data..."
                dists = np.load(fn_in)*nm2A
                steps = len(dists)
                dens = m/nm2A
                max_dist = np.max(dists)
                print "Loaded Data"
   
                if nd == 1: 
                    # Calculating the linear density

                    print "Calculating Linear Density"
                    nwaters = float(n)
                    vol_dens = nwaters/max_dist 
                    jacobian = 1.0

                elif nd == 2: 
                    # Calculating the area density
                    print "Calculating Area Density"
                    print "Code Segment is unfinished here! Edit at line 106"
                    exit()

                    #nwaters = float(n)
                    #vol_dens = nwaters/max_dist 
                    #jacobian = 1.0

                elif nd == 3: 
                    # Calculating the Volume density

                    print "Calculating Volume Density"
                    nwaters = float(n)
                    jacobian = 4.0*math.pi
   
                # Calculating Histogram Values 

                print "Calculating Histogram"
                bins = np.arange(0.0,max_dist,bin_interval)
                histVals,bin_edges = np.histogram(dists,bins=bins,density=False)
                histVals = np.array(histVals)
                histVals = histVals.astype(float)
    
                bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:]) 
                if nd == 2:
                    #TODO:
                    print "NEED to do 2 D RDF"
                    exit()
                elif nd == 3:
                    jacobian = np.square(bin_centers)*jacobian
                bin_width = bin_centers[1]-bin_centers[0]

                # => Controlling for Finite Wire Sampling <= #
                if periodic == False and nd == 1:

                    if fn_in == fn_Ostar_O or fn_in == fn_Hstar_O:
                        scale_val = float(n)-2.0
                        w_contribs = np.linspace(scale_val,1.0,len(bin_centers))/(scale_val/2.0)
                    elif fn_in == fn_Ostar_H:
                        scale_val = float(n)-2.0
                        w_contribs = np.linspace(scale_val,1.0,len(bin_centers))/(scale_val/4.0)
                    else:
                        scale_val = float(n)
                        w_contribs = np.linspace(scale_val,1.0,len(bin_centers))
                else: 
                #TODO: write finite w_contribs for 2 and 3 D systems. 
                    w_contribs = 1.0

                # Although normalizing by the density forces the RDF to converge to 1,
                # however, the relative peak heights are no longer meaningful when comparing
                # RDF's from systems with different var_1 

                if dens_norm == True:
                    #volume = (4.0/3.0)*math.pi*((max_dist/2.0)**3.0)
                    norm = steps*bin_width*vol_dens*jacobian*w_contribs
                else:
                    norm = steps*bin_width*jacobian*w_contribs

                histVals /= norm
    
                # => Testing that the RDF converges to 1 at long range <= #

                p = 0.5
                p2 = 0.25
                ind1 = int(len(histVals)*p)
                ind2 = int(len(histVals)*p2)

                mean = np.mean(histVals[ind2:ind1])
                ind3 = [indv for indv,val in enumerate(bin_centers[1:]) if val < 2.0 and val > 0.0]
                area = np.trapz(y=histVals[ind3],x=bin_centers[ind3])
                print "Peak area: %3.3f counts" % (area)
                print "Convergence: %3.4f" % (mean)
    
                # => Plotting the RDF <= #
    
                j = color_ind
                if ind == 0.0:
                    plt.plot(bin_centers[1:],histVals[1:],label=label1[ind],c=(0.0,1.0-j*ind,0.0+j*ind))
                elif ind == len(var_1)-1:
                    plt.plot(bin_centers[1:],histVals[1:],label=label1[ind],c=(0.0,1.0-j*ind,0.0+j*ind))
                else:
                    plt.plot(bin_centers[1:],histVals[1:],c=(0.0,1.0-j*ind,0.0+j*ind))
    
        	    os.chdir(ref_dir)

    os.chdir(ref_dir)
    hs.close()
   
    if ones_line == True: 
        x = np.arange(0.0,max_dist)
        y = np.ones(len(x))*3.3
        plt.plot(x,y,'--r',label=r'$\mathrm{Lit.\ \rho}$')

    plt.legend()

    if fn_in == fn_Hstar_O:
        plt.xlabel(r'$\mathrm{H^{*}-O\ distances\ \AA}$',fontsize=20)
    elif fn_in == fn_Ostar_O:
        plt.xlabel(r'$\mathrm{O^{*}-O\ distances\ \AA}$',fontsize=20)
    elif fn_in == fn_Ostar_H:
        plt.xlabel(r'$\mathrm{O^{*}-H\ distances\ \AA}$',fontsize=20)
    elif fn_in == fn_O_O:
        plt.xlabel(r'$\mathrm{O-O\ distances\ \AA}$',fontsize=20)
    elif fn_in == fn_O_H:
        plt.xlabel(r'$\mathrm{O-H\ distances\ \AA}$',fontsize=20)
    elif fn_in == fn_H_H:
        plt.xlabel(r'$\mathrm{H-H\ distances\ \AA}$',fontsize=20)

    if dens_norm == True:
        plt.ylabel(r'$\mathrm{g(r)}$',fontsize=25)
    else:
        plt.ylabel(r'$\mathrm{\rho\cdot g(r)}$',fontsize=25)

    if len(axis_range) != 1:
        plt.axis(axis_range)
    plt.savefig(fn_pdf)


def run_rdf():

    # Dump files for different types of RDF calculations 
   
     
    fn_Ostar_O = 'Ostar-O_dists' # Calc. Correct
    fn_Ostar_H = 'Ostar-H_dists'
    fn_Hstar_O = 'Hstar-O_dists' # Calc. Correct
    fn_O_O = 'O-O_dists' # Calc. Correct
    fn_O_H = 'O-H_dists'
    fn_H_H = 'H-H_dists'
    
    # => Densities <= #
    
    var_1 = [
                 3.0,
                 3.2,
                 3.4,
                 3.5,
                 3.6,
                 3.7,
                 3.8,
                 3.9,
                 4.0,
                 4.1,
                 4.2,
                 4.3,
                 4.4,
                 4.5,
                 ]

    waters = [100]

    #plot_rdf(
    #        fn_in=fn_Ostar_O,
    #        waters=waters,
    #        var_1=var_1,
    #        bin_interval = 0.04,
    #        ones_line = False,
    #        calculate_distances = False,
    #        ref_dir = '/Users/monikawilliams/Desktop/current_Research/var_dens/',
    #        color_ind = 0.06,
    #        axis_range = [0.0,20.0,0.0,10.0],
    #        dens_norm = False,
    #        )

    #plot_rdf(
    #        fn_in=fn_Hstar_O,
    #        waters=waters,
    #        var_1=var_1,
    #        bin_interval = 0.04,
    #        ones_line = False,
    #        calculate_distances = False,
    #        ref_dir = '/Users/monikawilliams/Desktop/current_Research/var_dens/',
    #        color_ind = 0.06,
    #        axis_range = [0.0,10.0,0.0,5.0],
    #        dens_norm = False,
    #        )

    #plot_rdf(
    #        fn_in=fn_O_O,
    #        waters=waters,
    #        var_1=var_1,
    #        bin_interval = 0.04,
    #        ones_line = False,
    #        calculate_distances = False,
    #        ref_dir = '/Users/monikawilliams/Desktop/current_Research/var_dens/',
    #        color_ind = 0.06,
    #        axis_range = [0.0,10.0,0.0,10.0],
    #        dens_norm = False,
    #        )

    #plot_rdf(
    #        fn_in=fn_Ostar_H,
    #        waters=waters,
    #        var_1=var_1,
    #        bin_interval = 0.04,
    #        ones_line = False,
    #        calculate_distances = False,
    #        ref_dir = '/Users/monikawilliams/Desktop/current_Research/var_dens/',
    #        color_ind = 0.06,
    #        axis_range = [0.0,10.0,0.0,10.0],
    #        dens_norm = False,
    #        )

#run_rdf()


    
    



