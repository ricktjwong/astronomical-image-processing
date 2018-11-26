##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.0
## Maintainer: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
##################################################

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

hdulist = fits.open("../data/mosaic.fits")
data = hdulist[0].data

"""
Manually mask the bright image
"""

def mask(data, x, y):
    x1, x2 = x[0], x[1]
    y1, y2 = y[0], y[1]
    for i in data[y1:y2]:
        i[x1:x2] = [0 for j in i[x1:x2]]
    return data

# Mask the circle
mask(data, [1226, 1626], [2990, 3430])
#plt.plot(data[3211])

# Mask the bleeding rectangles
mask(data, [1428, 1458], [3430, 4315])  # top lower rectangle
mask(data, [1434, 1445], [4315, 4611])  # top higher rectangle

mask(data, [1425, 1458], [1977, 2990])  # bottom higher rectangle
mask(data, [1430, 1444], [0, 1977])     # top lower rectangle

