##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.1
## Maintainers: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
## Description: Background fitting with a Gaussian
## distribution
##################################################

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
sys.path.append("../")
import utils.plot as pt

hdulist = fits.open("../data/fits/masked.fits")
data = hdulist[0].data
data = data.astype(np.float64)

centres = np.load("../data/detection/centres_12.npy")
galactic_intensities = np.load("../data/detection/galactic_intensities_12.npy")
background_intensities = np.load("../data/detection/background_intensities_12.npy")
img_data = np.load('../data/detection/final_img_12.npy')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
data = np.ma.masked_where(data < 10, data)
plt.imshow(data)
pt.mark_detected_objects(centres, ax)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
img_data = np.ma.masked_where(img_data < 10, img_data)
plt.imshow(img_data)
pt.mark_detected_objects(centres, ax)

z = 25.3
intensities = [i for i in galactic_intensities if i > 0]
m = -2.5 * np.log10(intensities) + z

def filter_m(m, threshold):
    return len([i for i in m if i < threshold])

x = np.linspace(m.min(), m.max())
y = [filter_m(m, i) for i in x]
plt.figure()
plt.plot(x, np.log10(y))