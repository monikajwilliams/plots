#!/usr/bin/env python
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import math
import pylab
import copy

##<=======================Set Values====================>
R = 8.314 #J/Kmol
R_kJ = R/1000
T = 300 
r_0 = 0.7
kappa = 13

##<=======================Energy Calculation====================>
r = np.arange(0.0, 2.0, 0.005)
N = math.sqrt((2*math.pi*R_kJ*T)/kappa)

r2 = copy.deepcopy(r)
for i in range(len(r2)):
   if r2[i] < r_0:
	r2[i] = r_0

Ei = 0.5*kappa*((r2-r_0)**2)
ni = np.exp(-Ei/(R_kJ*T))
ni_N = ni/N

a2i = -0.2281
a4i = 1.09
a6i = 0.2341
a8i = 0.3254

a2 = a2i/4.2
a4 = a4i/4.2
a6 = a6i/4.2
a8 = a8i/4.2

Ev = a2*(r**2) + a4*(r**4) + a6*(r**6) + a8*(r**8)
ni2 = np.exp(-Ev/(R_kJ*T))
ni_N2 = ni2/N

plt.figure()
#plt.plot(r,Ei, 'b')
#plt.plot(r,Ev, 'r')
plt.plot(r,ni_N, 'b')
plt.plot(r,ni_N2, 'r')
plt.axis([0.0, 2.0, 0.0, 3.0])
plt.show()


##<=======================Data====================>

def plot(filename):
	data = np.loadtxt(filename)
        rhos = data[50,1:]
	rhos *= 10
        plt.figure()
        hist, bins = np.histogram(rhos, density=True, bins=50)
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width, color='b')
        plt.plot(center,hist,'b')
	plt.plot(r,ni_N, 'r')
#	plt.axis([0.0, 1.0, 0.0, 10.0])
        plt.show()

#<=======================Main===================>
#plot = plot('rhos.txt')
