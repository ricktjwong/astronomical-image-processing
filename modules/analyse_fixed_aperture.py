##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.1
## Maintainers: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
## Description: Analyse fixed aperture detection
##################################################

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


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
r = [5, 10]
c = ['r', 'b']
plt.rcParams['figure.figsize'] = 7, 5
plt.figure()
for i in range(len(r)):
    centres = np.load("../data/detection_sweep/fixed_aperture/centres_"+ str(r[i]) + "radius.npy")
    galactic_intensities = np.load("../data/detection_sweep/fixed_aperture/galactic_intensities_" + str(r[i]) + "radius.npy")
    background_intensities = np.load("../data/detection_sweep/fixed_aperture/background_intensities_" + str(r[i]) + "radius.npy")
    img_data = np.load("../data/detection_sweep/fixed_aperture/final_img_" + str(r[i]) + "radius.npy")
    intensities = [i for i in galactic_intensities if i > 0]
    m = -2.5 * np.log10(intensities) + z
    x = np.linspace(m.min(), m.max())
    y = [filter_m(m, i) for i in x]
    x = x[1:]    
    y = y[1:]
    start_idx, end_idx = 11, 36
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
    plt.plot(x[start_idx:end_idx], y_fit, '-', linewidth=1.5, c=c[i])
    plt.scatter(x, np.log10(y), marker='.', s=15, c=c[i])
    print(len(centres))
    print("m" + str(i) + ": ")
    print(m)
    print("error: ")
    print(cov[0][0])
