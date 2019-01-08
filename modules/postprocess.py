# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 19:52:53 2019

@author: Daniel
"""

import numpy as np
import matplotlib.pyplot as plt

z = 25.3
intensities = np.load("../data/detection_sweep/galactic_intensities_1sigma.npy") 
intensities2 = np.load("../data/7Jan_2r_oldmask_3sig/galactic_intensities_3sigma.npy")
intensities = [i for i in intensities if i > 0]
m = -2.5 * np.log10(intensities) + z

#pad = np.full(1000, 18.)
#pad2 = np.full(500, 17.5)

def filter_m(m, threshold):
    return len([i for i in m if i < threshold])            

#x = np.linspace(m.min(), m.max())
#y = [filter_m(m, i) for i in x]
#plt.figure()
#plt.plot(x, np.log10(y))
#plt.show()   

bins = np.arange(10, 18.5, 0.5)
x = np.arange(11,22,0.5)
x2 = np.arange(11,18.5,0.5)
#x2 = np.arange(11,20.5,0.5)
y = [4,3,6,23,28,73,118,315,577,0, 1513,3050,5607,10532,18239,16890\
     ,29443,50206,80324,123213,181428,242743]
y2 = [4,3,6,23,28,73,118,315,577,1000, 1513,3050,5607,10532,18239]
#y2 = [4,3,6,23,28,73,118,315,577,0, 1513,3050,5607,10532,18239,16890\
#     ,29443,50206,80324]
y2 = np.array(y2)

total = np.sum(y2)
y2 = y2/total

#m2 = np.hstack((m,pad,pad2))

plt.figure()
plt.hist(m, bins, None, True)
plt.plot(x2,y2)
plt.show()

