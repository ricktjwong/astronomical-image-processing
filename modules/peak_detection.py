##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.0
## Maintainer: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
##################################################

from math import sqrt
from skimage.feature import blob_dog, blob_log, blob_doh
from astropy.io import fits
import matplotlib.pyplot as plt

hdulist = fits.open("../data/mosaic.fits")
img = hdulist[0].data

blobs_log = blob_log(img, max_sigma=30, num_sigma=10, threshold=.01)
# Compute radii in the 3rd column.
blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)

fig, ax = plt.subplots()
#coords_radius = []
#
for blob in blobs_log:
    ax.set_title("LOG")
    ax.imshow(img, interpolation='nearest')
    y, x, r = blob
#    coords_radius.append([x, y, r])
    c = plt.Circle((x, y), r, color='red', linewidth=2, fill=False)
    ax.add_patch(c)

plt.show()
