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

hdulist = fits.open("../data/mosaic.fits")
data = hdulist[0].data

for i in mk.CIRCLES:
    data = mk.mask_circle(data, i[0], i[1])

for i in mk.RECTS:
    data = mk.mask_rect(data, i)

plt.figure(1)
plt.imshow(data)

total_intensity = sp.circle_intensity(data, (2459, 3215), 14)
raw, r_h = sp.half_light_radius(data, (2457, 3211), 14, total_intensity)

data = data - 3421
y = [sp.surface_intensity(data, (2459, 3215), i) for i in range(1, 13)]
x = [i for i in range(1, len(y) + 1)]
plt.figure()
plt.plot(x, y)

I_e = sp.surface_intensity(data, (2459, 3215), int(r_h))

I = [sp.I(i, I_e, int(r_h), 0.60) for i in x]
plt.plot(x, I)

for j in range(1,10):
    I = [sp.I(i, I_e, int(r_h), j) for i in x]
    plt.plot(x, I)

plt.figure()
total_intensity = sp.circle_intensity(data, (518, 3539), 3)
raw, r_h = sp.half_light_radius(data, (518, 3539), 3, total_intensity)

y = [sp.surface_intensity(data, (518, 3539), i) for i in range(1, 4)]
x = [i for i in range(1, len(y) + 1)]
plt.figure()
plt.plot(x, y)

I_e = sp.surface_intensity(data, (518, 3539), int(r_h))

for j in range(1,10):
    I = [sp.I(i, I_e, int(r_h), j) for i in x]
    plt.plot(x, I)
