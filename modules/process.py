##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.0
## Maintainer: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
##################################################

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema

hdulist = fits.open("../data/mosaic.fits")
data = hdulist[0].data
print("data:")
print(data)
print (type(data))
print (len(data))
print (data[0])

def G(x, mu, sigma):
    print(mu)
    print(sigma)
    return (1/(2*np.pi)**0.5*sigma) * np.exp(-((x - mu)**2)/(2*sigma))

#def sigma(G):
#    mu = np.mean(G)
#    return (1/len(G) * sum([(i - mu)**2 for i in G])) ** 0.5

def get_brightest(nparray):
    max = 0
    for i in nparray:
        if i > max:
            max = i
    return max


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def mask(data, threshold=30000):
    c = len(data[0])
    r = len(data)
    for j in range(c):
        for i in range(r):
            if data[i][j] > threshold:
                data[i][j] = 0.
                

x = np.array([[1,1,1,1,1], [1,1,1,1,2]])

#for i in range(200):
#    print(get_brightest(data[i]))
#    plt.hist(i, 20)
    
data_flatten = data.flatten()

#plt.hist(data_flatten)

#one_row = data[3211]
#smoothened_row = smooth(one_row, 10)
#x_array = [i for i in range(len(data[0]))]
#plt.plot(x_array, one_row)
#plt.plot(x_array, smooth(one_row, 1))
#plt.plot(x_array, smooth(one_row, 5))
#plt.plot(x_array, smoothened_row)

#extremas = argrelextrema(smoothened_row, np.greater, order=10)
#
#plt.figure(2)
#gaussian = one_row[2459 - 50 : 2459 + 50]
#offset_gaussian = [i - min(gaussian) for i in gaussian]
#candidate_gaussian = smooth(offset_gaussian, 10)
#sigma = np.std(candidate_gaussian)
#x = np.linspace(2459 - 50, 2459 + 50, 100)
#gaussian_fit = G(x, 2459, sigma)
#
#area1 = np.trapz(candidate_gaussian)
#area2 = np.trapz(gaussian_fit)
#gaussian_fit2 = G(x, 2459, sigma) * area1/area2
#
#plt.plot(x, candidate_gaussian)
#plt.plot(x, gaussian_fit2)
#
#mask(data)
#x_array = [i for i in range(len(data[0]))]
#plt.plot(x_array, one_row)
