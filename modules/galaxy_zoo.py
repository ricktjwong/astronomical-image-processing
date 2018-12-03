#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 10:40:45 2018

@author: ricktjwong
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image

import sersic_profile as sp
import mask as mk
import process as ps
import peak_detection as pk

def rgb2gray(rgb):
    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b
    return gray

#hdulist = fits.open("../data/galaxyzoo1.fits")
hdulist = fits.open("../data/psField-001740-5-0022.fit")
data = hdulist[1].data
print(data.dtype)
print(len(data))
print("Uncertain galaxy:")
print(data[0])
print("Spiral galaxy:")
print(data[1])
print("Elliptical galaxy:")
print(data[3])

#img = mpimg.imread("../data/spiral2.jpeg")
#img = rgb2gray(img)
#
#plt.imshow(img)
#
#peaks = pk.find_peaks(img)
#
#total_intensity = sp.circle_intensity(img, (50, 51), 24)
#raw, r_h = sp.half_light_radius(img, (50, 51), 24, total_intensity)
#
#print("Half-radius:" + str(r_h))
#
#y = [sp.surface_intensity(img, (50, 51), i) for i in range(1, 10)]
#print(y)
#x = [i for i in range(1, len(y) + 1)]
#plt.figure()
#plt.plot(x, y, c='b')
#
#I_e = sp.surface_intensity(img, (50, 51), int(r_h))
#print(I_e)


