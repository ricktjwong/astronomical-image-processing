##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.0
## Maintainer: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
##################################################

from astropy.io import fits

hdulist = fits.open("../data/mosaic.fits")
data = hdulist[0].data

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

CIRCLES = [[(1436, 3214), 36], [(775, 3319), 12], [(558, 4095), 5],
           [(905, 2284), 12], [(2089, 1424), 7]]

RECTS = [[1418, 1454, 3236, 4315], [1428, 1444, 4315, 4610],
         [1418, 1454, 1977, 3200], [1430, 1444, 0, 1977],
         [900, 910, 2292, 2351], [900, 910, 2231, 2275],
         [767, 784, 3320, 3416], [767, 784, 3210, 3320]]
