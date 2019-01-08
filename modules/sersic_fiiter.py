# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 10:19:49 2018
@author: Daniel
"""

import sys
sys.path.append("../")
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import utils.plot as pt
import modules.peak_detection as pk
import sersic_profile as s_p

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 16
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.figsize'] = 5, 5
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'

hdulist2 = fits.open("../data/fits/mosaic.fits")
fulldata = hdulist2[0].data
fulldata = fulldata.astype(np.float64)
fulldata_subtracted = fulldata - 3450

hdulist = fits.open("../data/fits/masked.fits")
data = hdulist[0].data
data = data.astype(np.float64)

test = data[901:1101, 901:1101]
test = test - 3421

def sersic(Ie, Re, R, n):
    """ Ie is intensity at half light radius Re; Ie, Re are constant parameters. """
    b = 2*n**(-1/3)
    return Ie * np.exp(-b * ( (R/Re)**(1/n) - 1) )

def gauss(R, sigma):
    return 1/(np.sqrt(np.pi * 2) * sigma) * np.exp(-(R*R)/(2*sigma*sigma))

def moffat(r, R, background, I0, beta):
    return background + I0/(1 + (r*r)/(R*R) ) ** beta 

def moffat_fit(r, R, background, I0, I, guess):
    def moffat_beta(r, beta):
        return background + I0/(1 + (r*r)/(R*R) ) ** beta 
    coeff, var_matrix = curve_fit(moffat_beta, r, I, guess)
    return coeff, var_matrix

def sersic_fit(Ie, Re, I_data, R_data, guess):
    R_data, I_data = np.array(R_data), np.array(I_data)
    def sersic_n(R_data, n):
         b = 2*n**(-1/3)
         return Ie * np.exp(-b * ( (R_data/Re)**(1/n) - 1) )
    coeff, var_matrix = curve_fit(f = sersic_n, xdata = R_data, ydata = I_data, p0 = guess)
    return coeff, var_matrix

def gauss_fit(R,I,guess):
    normalisation = np.trapz(I,R)
    I_normalised = I / normalisation
    coeff, var_matrix = curve_fit(gauss, R, I_normalised, guess)
    return coeff, var_matrix, normalisation


sliced = fulldata[1380:1470, 2070:2110]
sliced = sliced - 3450
sliced2 = fulldata[2220:2350, 843:973] - 3450

surface_star = np.array([s_p.surface_intensity(sliced, (18, 46), r) for r in range(16)])
surface_star = surface_star / 10000
sigma_star, var, normalisation = gauss_fit([r for r in range(len(surface_star))], surface_star, 1)
beta_opt, moffat_var = moffat_fit([r for r in range(16)], 10, 0, 3.3550,  surface_star, 1)
print ("moffat params2: ", beta_opt, moffat_var)
moffat_star = [moffat(r, 10, 0, 3.3550, beta_opt) for r in range(16)]

plt.figure(5)
plt.imshow(sliced2)
plt.savefig("star.pdf", dpi=3000)
plt.show()

plt.figure(6)
#plt.title("star at centre(18,46)")
plt.scatter([r for r in range(16)], surface_star, c="black", marker='x')
plt.plot([r for r in range(16)], moffat_star, "r--")
plt.savefig("moffat.pdf", dpi=3000)
plt.show()
