#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 15:06:36 2018

@author: ricktjwong
"""

from math import sqrt
from skimage.feature import blob_dog, blob_log, blob_doh
from astropy.io import fits
import matplotlib.pyplot as plt
import mask as mk

hdulist = fits.open("../data/mosaic.fits")
data = hdulist[0].data

data = data - 3421

for i in mk.CIRCLES:
    data = mk.mask_circle(data, i[0], i[1])

for i in mk.RECTS:
    data = mk.mask_rect(data, i)

image = data
image_gray = data

blobs_log = blob_log(image_gray, max_sigma=30, num_sigma=10, threshold=.01)

# Compute radii in the 3rd column.
blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)

blobs_dog = blob_dog(image_gray, max_sigma=30, threshold=.001)
blobs_dog[:, 2] = blobs_dog[:, 2] * sqrt(2)

blobs_doh = blob_doh(image_gray, max_sigma=30, threshold=.01)

blobs_list = [blobs_log, blobs_dog, blobs_doh]
colors = ['yellow', 'lime', 'red']
titles = ['Laplacian of Gaussian', 'Difference of Gaussian',
          'Determinant of Hessian']
sequence = zip(blobs_list, colors, titles)

fig, axes = plt.subplots(1, 3, figsize=(9, 3), sharex=True, sharey=True)
ax = axes.ravel()

for idx, (blobs, color, title) in enumerate(sequence):
    ax[idx].set_title(title)
    ax[idx].imshow(image, interpolation='nearest')
    for blob in blobs:
        y, x, r = blob
        c = plt.Circle((x, y), r, color=color, linewidth=2, fill=False)
        ax[idx].add_patch(c)
    ax[idx].set_axis_off()

plt.tight_layout()
plt.show()