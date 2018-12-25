##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.0
## Maintainer: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
##################################################

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

hdulist = fits.open("../data/masked.fits")
data = hdulist[0].data


# Define model function to be used to fit to the data above:
def gauss(x,*p):
    mu , sigma = p
    return (1/(sigma*np.sqrt(2*np.pi)) * np.exp(-(x-mu)**2/(2.*sigma**2)))

def removeHighValues(data):
    new_data = np.zeros(data.shape)
    for i in range(len(data)):
        if data[i] > 5000:
            new_data[i] = 0
        else:
            new_data[i] = data[i]
    return new_data            

plt.figure(1, figsize = (6,8))
plt.imshow(data)
plt.show()
                 
noise_data = data.flatten()
bins = np.linspace(3200, 3700, 100)

plt.figure()
hist = plt.hist(noise_data, bins, None, True)
plt.show()
counts, intensity_bins = hist[0], hist[1]
bin_centre = (intensity_bins[:-1] + intensity_bins[1:])/2

p0 = [3421, 10]
coeff, var_matrix = curve_fit(gauss, bin_centre, counts, p0=p0)
