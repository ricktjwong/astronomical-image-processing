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

R_range = [r for r in range(2,21)]
R_range2 = [r for r in range(2,16)]

total_I, total_pix =  s_p.circle_intensity(test, (65,122), 20)
Re, radius_corrected = s_p.half_light_radius(test, (65,122), 20, total_I) # half light radius
surface_I = [s_p.surface_intensity(test, (65, 122), r) for r in range(21)]
Ie = surface_I[Re]

total_I, total_pix =  s_p.circle_intensity(test, (133,26), 15)
Re2, radius_corrected2 = s_p.half_light_radius(test, (133,26), 15, total_I) # half light radius
surface_I2 = [s_p.surface_intensity(test, (133, 26), r) for r in range(16)]
Ie2 = surface_I[Re2]

n_opt, var_matrix = sersic_fit(Ie, Re, surface_I[2:], R_range, 1)
print (n_opt, var_matrix)
n_opt2, var_matrix2 = sersic_fit(Ie2, Re2, surface_I2[2:], R_range2, 1)
print (n_opt2, var_matrix2)

sigma_opt, var, normalisation = gauss_fit([r for r in range(len(surface_I2))], surface_I2, 1)
print (sigma_opt, var, normalisation)

gauss_fitted = [normalisation * gauss(r, sigma_opt) for r in range(16)]

sersic_fitted = [sersic(Ie, Re, r, n_opt) for r in range(21)]

plt.figure(1)
plt.imshow(test)
plt.show()

surface = [s_p.surface_intensity(test, (65, 122), r) for r in range(21)]
sigma_star, var, normalisation = gauss_fit([r for r in range(len(surface))], surface, 1)
#print ("gauss parameters: ", (sigma_star, var, normalisation))
gauss_star = [normalisation * gauss(r, sigma_star) for r in range(21)]
beta_opt, moffat_var = moffat_fit([r for r in range(21)], 15, 0, 500,  surface, 1)
print ("moffat params0: ", beta_opt, moffat_var)
moffat_star = [moffat(r ,15, 0, 500, beta_opt) for r in range(21)]

plt.figure(2)
plt.title("centre at (65,122)")
plt.plot(range(len(surface_I)), surface_I, "b.")
plt.plot(R_range, sersic_fitted[2:], "yx")
plt.plot(range(len(moffat_star)), moffat_star, "go")
plt.show()

sersic_fitted = [sersic(Ie, Re, r, n_opt) for r in range(16)]
surface_star = [s_p.surface_intensity(test, (133, 26), r) for r in range(16)]
beta_opt, moffat_var = moffat_fit([r for r in range(16)], 10, 0, 1020,  surface_star, 1)
print ("moffat params1: ", beta_opt, moffat_var)
moffat_star = [moffat(r ,15, 0, 1020, beta_opt) for r in range(16)]

plt.figure(3)
plt.title("centre at (133,26)")
plt.plot(range(len(surface_I2)), surface_I2, "b.")
plt.plot(R_range2, sersic_fitted[2:], "yx")
plt.plot(range(len(gauss_fitted)), gauss_fitted, "rx")
plt.plot(range(len(moffat_star)), moffat_star, "gx")
plt.show()

sliced = fulldata[1380:1470, 2070:2110]
sliced = sliced - 3450
sliced2 = fulldata [2220:2350, 875:925] - 3450

surface_star = [s_p.surface_intensity(sliced, (18, 46), r) for r in range(16)]
sigma_star, var, normalisation = gauss_fit([r for r in range(len(surface_star))], surface_star, 1)
#print ("gauss parameters: ", (sigma_star, var, normalisation))
gauss_star = [normalisation * gauss(r, sigma_star) for r in range(16)]
beta_opt, moffat_var = moffat_fit([r for r in range(16)], 10, 0, 33550,  surface_star, 1)
print ("moffat params2: ", beta_opt, moffat_var)
moffat_star = [moffat(r ,10, 0, 33550, beta_opt) for r in range(16)]

plt.figure(4)
plt.imshow(fulldata)
plt.show()

plt.figure(5)
plt.imshow(sliced2)
plt.show()

plt.figure(6)
plt.title("star at centre(18,46)")
plt.plot([r for r in range(16)], surface_star, "b.")
plt.plot([r for r in range(16)], gauss_star, "yx")
plt.plot([r for r in range(16)], moffat_star, "rx")
plt.show()

surface_star2 = [s_p.surface_intensity(sliced2, (29, 66), r) for r in range(21)]
sigma_star, var, normalisation = gauss_fit([r for r in range(len(surface_star2))], surface_star2, 1)
print ("gauss parameters: ", (sigma_star, var, normalisation))
gauss_star2 = [normalisation * gauss(r, sigma_star) for r in range(21)]
beta_opt, moffat_var = moffat_fit([r for r in range(21)], 15, 0, 36800,  surface_star2, 1)
print ("moffat params3: ", beta_opt, moffat_var)
moffat_star2 = [moffat(r ,15, 0, 36800, beta_opt) for r in range(21)]

plt.figure(7)
plt.title("star at centre(29,66)")
plt.plot([r for r in range(21)], surface_star2, "b.")
plt.plot([r for r in range(21)], gauss_star2, "yx")
plt.plot([r for r in range(21)], moffat_star2, "rx")
plt.show()

surface_star3 = [s_p.surface_intensity(fulldata_subtracted, (966, 1654), r) for r in range(11)]
sigma_star, var, normalisation = gauss_fit([r for r in range(len(surface_star3))], surface_star3, 1)
print ("gauss parameters: ", (sigma_star, var, normalisation))
gauss_star3 = [normalisation * gauss(r, sigma_star) for r in range(11)]
beta_opt, moffat_var = moffat_fit([r for r in range(11)], 10, 0, 32850,  surface_star3, 1)
print ("moffat params4: ", beta_opt, moffat_var)
moffat_star3 = [moffat(r ,10, 0, 32850, beta_opt) for r in range(11)]

total_I, total_pix =  s_p.circle_intensity(fulldata_subtracted, (966, 1654), 10)
Re, radius_corrected = s_p.half_light_radius(fulldata_subtracted, (966, 1654), 10, total_I)
Ie = surface_star3[Re]
n_opt, var_matrix = sersic_fit(Ie, Re, surface_star3[2:], [r for r in range(2,11)], 1)
print ("sersic params: ", n_opt, var_matrix)
sersic_fitted = [sersic(Ie, Re, r, n_opt) for r in range(0,11)]

plt.figure(8)
plt.title("centre(966,1654)")
plt.plot([r for r in range(11)], surface_star3, "b.")
plt.plot([r for r in range(0,11)], sersic_fitted, "yx")
plt.plot([r for r in range(11)], gauss_star3, "rx")
plt.plot([r for r in range(11)], moffat_star3, "gx")
plt.show()

surface_star3 = [s_p.surface_intensity(fulldata_subtracted, (775, 3320), r) for r in range(21)]
sigma_star, var, normalisation = gauss_fit([r for r in range(len(surface_star3))], surface_star3, 1)
#print ("gauss parameters: ", (sigma_star, var, normalisation))
gauss_star3 = [normalisation * gauss(r, sigma_star) for r in range(21)]
beta_opt, moffat_var = moffat_fit([r for r in range(21)], 15, 0, 36550,  surface_star3, 1)
print ("moffat params5: ", beta_opt, moffat_var)
moffat_star3 = [moffat(r, 15, 0, 36550, beta_opt) for r in range(21)]

plt.figure(9)
plt.title("star at centre(775,3320)")
plt.plot([r for r in range(21)], surface_star3, "b.")
plt.plot([r for r in range(21)], gauss_star3, "yx")
plt.plot([r for r in range(21)], moffat_star3, "rx")
plt.show()
