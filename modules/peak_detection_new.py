#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 14:27:33 2018

@author: ricktjwong
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import math

import sersic_profile as sp
import mask as mk
import process as ps
import peak_detection as pk

z = 25.3

hdulist = fits.open("../data/mosaic.fits")
data = hdulist[0].data
data = data.astype(np.float64)

for i in mk.CIRCLES:
    data = mk.mask_circle(data, i[0], i[1])

for i in mk.RECTS:
    data = mk.mask_rect(data, i)
    
#data = data[101:201, 101:201]
rows = data.shape[0]
cols = data.shape[1]

max_flattened = np.argmax(data)
max_r = int(max_flattened / cols)
max_c = max_flattened % cols

print(max_r, max_c)
print(data[max_r, max_c])

plt.figure(1)
plt.imshow(data)

r = 5
threshold = 3485

box_x1 = max_c - r
box_x2 = max_c + r
box_y1 = max_r - r    
box_y2 = max_r + r
final = []
n = 0
galactic_intensities = []
background_intensities = []
#while (max_flattened != 0):
while (n < 10000):
    galaxy_intensity = []
    background_intensity = []
    box_x1 = max_c - r
    box_x2 = max_c + r
    box_y1 = max_r - r    
    box_y2 = max_r + r
    print(box_x1, box_x2, box_y1, box_y2)
    if box_x1 < 1 or box_x2 > cols-1 or box_y1 < 1 or box_y2 > rows-1:
        data[max_r][max_c] = 0
        max_flattened = np.argmax(data)
        max_r = int(max_flattened / cols)
        max_c = max_flattened % cols        
        n += 1
    else:
        for i in range(box_y1, box_y2 + 1):
            for j in range(box_x1, box_x2 + 1):
                if data[i][j] < threshold:
                    background_intensity.append(data[i][j])
                else:
                    galaxy_intensity.append(data[i][j])
                    data[i][j] = 0
                    
#        if len(background_intensity) < 5:
#            r += 2
            # repeat
        galactic_intensities += [sum(galaxy_intensity)]
        background_intensities += [sum(background_intensity)]
        max_flattened = np.argmax(data)
        max_r = int(max_flattened / cols)
        max_c = max_flattened % cols
        n += 1

plt.figure(2)
plt.imshow(data)

