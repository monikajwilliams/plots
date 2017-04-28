#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

def plot_PMF(axis = 2,
            fn_in = "Ostar",
            fn_PMF ="PMF.pdf",
            fn_PDF ="PDF.pdf",
            num_bins = 25.0,
            T = 300.0,
            ):


    # => Loading data <= #    

    all_positions = np.load(fn_in)
    positions = all_positions[:,axis] - min(all_positions[:,axis])
    max_pos =  max(positions)
    print "Loaded Data"
   
    # => Plotting probability distribution <= # 

    plt.clf() 
    hist, bins = np.histogram(positions,num_bins, density=True)
    width = 0.7 * (bins[1] - bins[0])
    center = 0.5 * (bins[:-1] + bins[1:]) 
    
    plt.bar(center, hist, align='center', width=width)
    plt.xlabel(r"$\mathrm{Position\ \AA}$")
    plt.ylabel(r"$\mathrm{Probability}$")
    plt.savefig(fn_PDF)
    
    # => Calculating potential of mean force <= # 

    E = -np.log(hist)
    num_E = len(E)
    E_0 = min(E)
    E -= E_0
    x = np.linspace(0.0,max_pos,num=num_E)
    y = np.linspace(0.0,10.0,num=num_E)
    wire_end = np.ones(len(y))*max_pos
    
    # => Plotting potential of mean force <= # 

    plt.clf()
    plt.plot(x,E)
    plt.plot(wire_end,y,label=r"$\mathrm{wire\ length\ (\AA)}$")
    plt.xlabel(r"$\mathrm{Position\ (\AA)}$",fontsize=18)
    plt.ylabel(r"$\mathrm{Free\ Energy\ (k_{B}T)}$",fontsize=18)
    plt.legend(fontsize=14)
    plt.axis([0.0,max(x) + 0.1*max(x),0.0,10.0])
    plt.savefig(fn_PMF)

    return


def plot_Z(
           # Traj. timestep (fs)
           timestep = 50.0,
           # Time to drop from start of traj. (ps)
           drop_ps = 10.0,
           # Filename to load
           fn_in = "Ostar",
           # Filename to output
           fn_out = "Zs.png"
           ):
    
    # => Unit conversions <= #

    fs2ps = 0.001 
    ps2ns = 0.001

    # In timestep/frame units
    drop_traj = drop_ps/(timestep*fs2ps) 
    # Index for new traj. start
    ind_drop_traj = int(drop_traj) 
    
    # => Trajectory Relevant Data <= #

    Ostars = np.load(fn_in)
    OstarZs = Ostars[ind_drop_traj:,2]
    OstarZs -= min(OstarZs) 
    max_pos = max(OstarZs)
    print "Loaded Data"
    
    # Setting time
    timesL = np.arange(len(OstarZs))
    times = np.array(timesL*timestep*fs2ps)
    wire_end = np.ones(len(times))*max_pos
  
    plt.clf() 
    plt.plot(times,OstarZs,c='b')
    plt.plot(times,wire_end,label=r"$\mathrm{wire\ length\ (\AA)}$")
    plt.xlabel(r'$\mathrm{t\ (ps)}$',fontsize=18)
    plt.ylabel(r'$\mathrm{Z\ position\ O^*\ (\AA)}$',fontsize=18)
    plt.legend()
    plt.axis([0.0,max(times)+10.0,0.0,max_pos+10.0])
    plt.savefig(fn_out) 

    return
