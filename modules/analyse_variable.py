##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.1
## Maintainers: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
## Description: Background fitting with a Gaussian
## distribution
##################################################

import sys
sys.path.append("../")
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import utils.plot as pt


def filter_m(m, threshold):
    return len([i for i in m if i < threshold])


def polyfit(x, y, degree, weight):
    results={}
    coeffs,V=np.polyfit(x,y,degree,w=weight,cov=True)
    results['polynomial']=coeffs.tolist()                           # polynomial coefficients
    results['cov']=V.tolist()   
    p=np.poly1d(coeffs)                                             # r-squared computation starts here
    yHat=p(x)                                                       # fit values, and mean
    yBar=np.mean(y)  
    ssReg=np.sum((yHat-yBar)**2)   
    ssRot=np.sum((y - yBar)**2)
    results['determination']=ssReg/ssRot
    return results


hdulist = fits.open("../data/fits/masked.fits")
data = hdulist[0].data
data = data.astype(np.float64)
mu, sigma = 3418.5636925, 12.58387389
z = 25.3

#plt.hist(centres[:,2], bins=max(centres[:,2]) - min(centres[:,2]) + 1)
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#data = np.ma.masked_where(data < 10, data)
#plt.imshow(data)
#pt.mark_detected_objects(centres, ax)
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#img_data = np.ma.masked_where(img_data < 10, img_data)
#plt.imshow(img_data)
#pt.mark_detected_objects(centres, ax)

plt.figure(1)
for i in range(3, 4):
    centres = np.load("../data/detection_sweep/variable_aperture_r5/centres_" + str(i) + "sigma.npy")
    galactic_intensities = np.load("../data/detection_sweep/variable_aperture_r5/galactic_intensities_" + str(i) + "sigma.npy")
    background_intensities = np.load("../data/detection_sweep/variable_aperture_r5/background_intensities_" + str(i) + "sigma.npy")
    img_data = np.load("../data/detection_sweep/variable_aperture_r5/final_img_" + str(i) + "sigma.npy")
    intensities = [i for i in galactic_intensities if i > 0]
    m = -2.5 * np.log10(intensities) + z
    x = np.linspace(m.min(), m.max())
    y = [filter_m(m, i) for i in x]
    x = x[1:]    
    y = y[1:]
    start_idx, end_idx = 15, 31
    w = np.array([1 / y[i]**0.5 for i in range(len(x))])
    wx = [0.02 for i in range(len(w))]
    data = polyfit(x[start_idx:end_idx], np.log10(y[start_idx:end_idx]),
                   1, 1/w[start_idx:end_idx])
    coeffs = data['polynomial']
    cov = np.array(data['cov']) ** 0.5
    m = coeffs[0]
    b = coeffs[1]
    y_fit = [m * i + b for i in x[start_idx:end_idx]]
    plt.errorbar(x, np.log10(y),
                 yerr=w, linestyle="None", color='grey', linewidth=0.5)
    plt.errorbar(x, np.log10(y),
                 xerr=wx, linestyle="None", color='grey', linewidth=0.5)
    plt.plot(x[start_idx:end_idx], y_fit, '-', linewidth=1.5, c='g')
    plt.scatter(x, np.log10(y), marker='.', s=15, c='g')
    print(len(centres))
    print("m1: ")
    print(m)
    print("error: ")
    print(cov)
