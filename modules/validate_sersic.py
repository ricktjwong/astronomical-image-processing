##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.0
## Maintainers: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
## Description: Validate that galaxies have an
## intensity distribution well modelled by the
## Sersic profile using an actual galaxy from DR12
##################################################

from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import sersic_profile as sp
from profiling import sersic_fit

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 16
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.figsize'] = 5, 5
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'

hdulist = fits.open("../data/fits/dr12.fits")
data = hdulist[0].data
plt.figure()
plt.imshow(data, cmap='rainbow')
plt.scatter(974, 396, c="r", s=1)

"""
Galaxy from SDSS DR12
"""

w = WCS("../data/fits/dr12.fits")
px, py = w.all_world2pix(179.80335889, -0.523868, 1)
print(px, py)

x1, x2 = 974, 396
plt.figure()
plt.imshow(data[386:407, 965:986])
r = 10
total_intensity = sp.circle_intensity(data, (x1, x2), r)[0]
raw, r_h = sp.half_light_radius(data, (x1, x2), r, total_intensity)

data = data + 0.1315918
y = [sp.surface_intensity(data, (x1, x2), i) for i in range(0, r)]
print(y)
x = [i for i in range(len(y))]
plt.figure()
plt.scatter(x, y, c='black', marker='x')

I_e = sp.surface_intensity(data, (x1, x2), int(r_h))
print(I_e)
print(r_h)

for j in range(1, 4):
    I = [sp.I(i, I_e, int(r_h), j) for i in x]
#    plt.plot(x, I, "--")
    
#total_I, total_pix = sp.circle_intensity(data, (x1, x2), r)
Re, Re_corrected = sp.half_light_radius(data, (x1, x2), r, total_intensity)           # Half light radius
I_surface = [sp.surface_intensity(data, (x1, x2), i) for i in range(r+1)]    # Surface intensities (ydata)
Ie = I_surface[Re]                                                                  # Half light Intensity
n_opt, sersic_cov = sersic_fit(Ie, Re, I_surface, range(0, r+1), 1)              # Obtain fitted sersic index n

I = [sp.I(i, I_e, int(r_h), n_opt) for i in x]
plt.plot(x, I, "--", c='r')
plt.savefig("validate_sersic.pdf", dpi=3000)

"""
Star
"""

#px, py = w.all_world2pix(179.71819, -0.45674, 1)
#print(px, py)
#
#x1, x2 = 1645, 398
#r = 12
#total_intensity = sp.circle_intensity(data, (x1, x2), r)[0]
#raw, r_h = sp.half_light_radius(data, (x1, x2), r, total_intensity)
#
#y = [sp.surface_intensity(data, (x1, x2), i) for i in range(0, r)]
#print(y)
#x = [i for i in range(0, len(y))]
#plt.figure()
#plt.plot(x, y, c='b')
#
#I_e = sp.surface_intensity(data, (x1, x2), int(r_h))
#print(I_e)
#print(r_h)
#
#for j in range(1,4):
#    I = [sp.I(i, I_e, int(r_h), j) for i in x]
#    plt.plot(x, I, "--")
#
#x1, x2 = 1645, 398
#surf_i = []
#for i in range(12):
#    surf_i.append(data[x2][x1 + i])
#x = [i for i in range(len(surf_i))]
#plt.figure()
#plt.plot(x, surf_i)
#
#
#x1, x2 = 973, 397
#surf_i = []
#for i in range(12):
#    surf_i.append(data[x2][x1 + i])
#x = [i for i in range(len(surf_i))]
#plt.figure()
#plt.plot(x, surf_i)
