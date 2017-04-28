#!/usr/bin/env python
import matplotlib 
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import pandas as pd
import math
import pylab


def plot(filename,figurename,col):
	data = np.loadtxt(filename)
        x = data[:,1]
        fig = plt.figure()
        plt.plot(x,data[:,col])
	fig.savefig(figurename)

