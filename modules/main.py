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

z = 25.3

hdulist = fits.open("../data/masked.fits")
data = hdulist[0].data
data = data.astype(np.float64)

# Sersic profile fitting on an existing galaxy centred at 2459, 3215
plt.imshow(data[3205:3225, 2450:2470])
total_intensity = sp.circle_intensity(data, (2459, 3215), 14)[0]
raw, r_h = sp.half_light_radius(data, (2457, 3211), 14, total_intensity)

y = [sp.surface_intensity(data, (2459, 3215), i) for i in range(1, 13)]
print(y)
x = [i for i in range(1, len(y) + 1)]
plt.figure()
plt.plot(x, y, c='b')

I_e = sp.surface_intensity(data, (2459, 3215), int(r_h))
print(I_e)
print(r_h)

I = [sp.I(i, I_e, int(r_h), 0.60) for i in x]
plt.plot(x, I, "--")

for j in range(1,10):
    I = [sp.I(i, I_e, int(r_h), j) for i in x]
    plt.plot(x, I, "--")


#rects = [[20, 30, 20, 4611], [2560, 2570, 20, 4611], [20, 2570, 20, 30],
#         [20, 2570, 4601, 4611], [1198, 1654, 426, 490], 
#         [1078, 1686, 314, 362], [1390, 1486, 218, 260],
#         [1310, 1526, 114, 170], [1134, 1654, 20, 50],
#         [1418, 1454, 3236, 4315], [1418, 1454, 4315, 4610],
#         [1418, 1454, 1977, 3200], [1420, 1454, 20, 1977],
#         [900, 910, 2292, 2351], [900, 910, 2231, 2275],
#         [767, 784, 3320, 3416], [767, 784, 3210, 3320], 
#         [960, 990, 2700, 2840]]
#remove_indices = []
#for i in range(len(peaks)):
#    for j in rects:
#        if peaks[i][0] <= j[3] + 20 and peaks[i][0] >= j[2] - 20 and \
#        peaks[i][1] <= j[1] + 20 and peaks[i][1] >= j[0] - 20:
#            remove_indices.append(i)
#peaks = np.delete(peaks, remove_indices, 0)
#print(len(peaks))
#y = peaks[:,0]
#x = peaks[:,1]
#plt.scatter(x, y, s=2, c='r')
#
#intensities = []
#for i in peaks:
#    y, x, r = i
#    r = math.ceil(r)
#    new_centre = 2 * r
#    x_min, x_max = int(x - new_centre), int(x + new_centre)
#    y_min, y_max = int(y - new_centre), int(y + new_centre)
#    
#    if (x_min < 0 or x_max > 2569 or y_min < 0 or y_max > 4610):
#        continue
#    else:
#        # Mask the square containing the galactic object
#        square = data[y_min:y_max+1, x_min:x_max+1].copy()
#        mask = np.ones((new_centre*2 + 1, new_centre*2 + 1))
#        mask[r:2*r+2, r:2*r+2] = 0
#        masked = np.ma.array(square, mask = mask)
#        masked = np.ma.masked_where(masked > 4000, masked)
#        background = np.mean(masked)
#        for i in range(len(square)):
#            for j in range(len(square[0])):
#                if square[i][j] >= background:
#                    square[i][j] -= background
#        total_intensity, total = sp.circle_intensity(square, (new_centre, new_centre), r)
#        intensities.append(total_intensity)
#        
##intensities = np.loadtxt("intensities.txt")
#intensities = [i for i in intensities if i > 0]
#m = -2.5 * np.log10(intensities) + z
#
#def filter_m(m, threshold):
#    return len([i for i in m if i < threshold])
#
#x = np.linspace(m.min(), m.max())
#y = [filter_m(m, i) for i in x]
#plt.figure()
#plt.plot(x, np.log10(y))
