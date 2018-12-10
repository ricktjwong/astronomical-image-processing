##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.0
## Maintainer: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
##################################################

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

import sersic_profile as sp
import mask as mk
import process as ps
import peak_detection as pk

hdulist = fits.open("../data/1000-1101.fits")
data = hdulist[0].data
data = data.astype(np.float64)

plt.imshow(data)

peaks = pk.find_peaks(data)

rm_indices = []
for i in range(len(peaks)):
    if peaks[i][0] == 100. or peaks[i][0] == 0. or peaks[i][1] == 100. \
    or peaks[i][1] == 0:
        rm_indices.append(i)
peaks = np.delete(peaks, rm_indices, 0)

y = peaks[:,0]
x = peaks[:,1]
plt.scatter(x, y, s=2, c='r')

print(peaks)
print(len(peaks))
