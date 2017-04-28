#!/usr/bin/env python
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import math
import pylab

natoms = 32

#<=======================Data====================>

def plot(filename):
	data = np.loadtxt(filename)
        x = data[:,0]
        plt.figure()
        for i in range(2):
            plt.plot(x,data[:,i+1])
        plt.show()

#<=======================Main===================>
#plot = plot('rhos.txt')
