##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.1
## Maintainers: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
## Description: Background fitting with a Gaussian
## distribution
##################################################

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 16
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.figsize'] = 8, 5
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'

hdulist = fits.open("../data/fits/masked.fits")
data = hdulist[0].data


def gauss(x, *p):
    """
    Gaussian function to fit to data
    """
    mu, sigma = p
    return (1 / (sigma * np.sqrt(2 * np.pi)) *
            np.exp(-(x - mu) ** 2 / (2. * sigma ** 2)))

#new_data = list(filter(lambda x: x < 5000, data))

def removeHighValues(data):
    new_data = np.zeros(data.shape)
    for i in range(len(data)):
        if data[i] > 5000:
            new_data[i] = 0
        else:
            new_data[i] = data[i]
    return new_data


noise_data = data.flatten()
bins = np.linspace(3200, 3700, 501)

plt.figure()
hist = plt.hist(noise_data, bins, None, True, color='b', edgecolor='black')
counts, intensity_bins = hist[0], hist[1]
bin_centre = (intensity_bins[:-1] + intensity_bins[1:]) / 2

p0 = [3421, 10]
coeff, var_matrix = curve_fit(gauss, bin_centre, counts, p0=p0)
plt.plot(bins, gauss(bins, *coeff), 'r--')

new_counts = np.zeros(435)
new_counts[0:218] = counts[0:218]
new_counts[218:] = counts[0:217][::-1]

new_intensity_bins = intensity_bins[0:436]
counts = new_counts
intensity_bins = new_intensity_bins
bin_centre = (intensity_bins[:-1] + intensity_bins[1:]) / 2

p0 = [3421, 10]
coeff, var_matrix = curve_fit(gauss, bin_centre, counts, p0=p0)
plt.plot(bins, gauss(bins, *coeff), 'g--')

plt.xlim([3350, 3550])
print(coeff, var_matrix)
print("Error: ")
print(np.sqrt(np.diag(var_matrix)))
