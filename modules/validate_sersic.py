#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 10:40:45 2018

@author: ricktjwong
"""

from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt

import sersic_profile as sp
import peak_detection as pk

hdulist = fits.open("../data/dr12.fits")
w = WCS("../data/dr12.fits")
px, py = w.all_world2pix(179.80335889, -0.523868, 1)
print(px, py)

data = hdulist[0].data

plt.figure(1)
plt.imshow(data)

x1, x2 = int(px), int(py)
x1, x2 = 973, 397
r = 10
total_intensity = sp.circle_intensity(data, (x1, x2), r)
raw, r_h = sp.half_light_radius(data, (x1, x2), r, total_intensity)

data = data + 0.1315918
y = [sp.surface_intensity(data, (x1, x2), i) for i in range(1, r)]
print(y)
x = [i for i in range(1, len(y) + 1)]
plt.figure()
plt.plot(x, y, c='b')

I_e = sp.surface_intensity(data, (x1, x2), int(r_h))
print(I_e)
print(r_h)

for j in range(1,4):
    I = [sp.I(i, I_e, int(r_h), j) for i in x]
    plt.plot(x, I, "--")
