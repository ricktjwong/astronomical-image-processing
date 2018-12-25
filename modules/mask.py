##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.0
## Maintainer: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
##################################################

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

hdulist = fits.open("../data/mosaic.fits")
data = hdulist[0].data
data = data.astype(np.float64)

"""
Manually mask the bright image
"""

def mask_rect(data, coords):
    x1, x2, y1, y2 = coords[0], coords[1], coords[2], coords[3]
    for i in data[y1:y2]:
        i[x1:x2] = [0 for j in i[x1:x2]]
    return data


def mask_circle(data, centre, radius):
    max_y, min_y = centre[1] + radius, centre[1] - radius
    max_x, min_x = centre[0] + radius, centre[0] - radius    
    for i in range(min_y, max_y + 1):
        for j in range(min_x, max_x + 1):
            if ((i - centre[1]) ** 2  + (j - centre[0])**2) <= (radius * radius):
                data[i][j] = 0                
    return data


def save_to_fits(data):
    hdu = fits.PrimaryHDU(data)
    hdul = fits.HDUList([hdu])
    hdul.writeto('masked.fits')


CIRCLES = [[(1436, 3214), 170], [(775, 3319), 30], [(558, 4095), 7],
           [(905, 2284), 30], [(2089, 1424), 7], [(2466, 3412), 14],
           [(2133, 3758), 20], [(972, 2772), 30], [(2131, 2308), 16],
           [(1455, 4032), 20]]

RECTS = [[0, 100, 0, 4611], [2470, 2570, 0, 4611], [0, 2570, 0, 100],
         [0, 2570, 4511, 4611], [1198, 1654, 426, 490], 
         [1078, 1686, 314, 362], [1390, 1486, 218, 260],
         [1310, 1545, 114, 170], [1134, 1654, 0, 50],
         [1418, 1454, 3236, 4315], [1418, 1454, 4315, 4610],
         [1418, 1454, 1977, 3200], [1420, 1454, 0, 1977],
         [900, 910, 2220, 2358], [767, 784, 3200, 3416],
         [960, 990, 2700, 2840], [2125, 2140, 2284, 2334],
         [2460, 2470, 3382, 3440], [2127, 2140, 3705, 3805],
         [554, 563, 4080, 4118]]

for i in CIRCLES:
    data = mask_circle(data, i[0], i[1])

for i in RECTS:
    data = mask_rect(data, i)

masked = np.ma.masked_where(data < 10, data)
plt.figure()
plt.imshow(masked)
#save_to_fits(data)
