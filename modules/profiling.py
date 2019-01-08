# -*- coding: utf-8 -*-
"""
Created on Sun Dec 30 10:15:28 2018

@author: Daniel
"""

""" Importing Modules """
import sys
sys.path.append("../")
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sersic_profile as s_p

""" Function Definitions """
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
         b = 2 * n ** (-1/3)
         return Ie * np.exp(-b * ( (R_data/Re)**(1/n) - 1) )
    coeff, var_matrix = curve_fit(f = sersic_n, xdata = R_data, ydata = I_data, p0 = guess)
    return coeff, var_matrix

def gauss_fit(R,I,guess):
    normalisation = np.trapz(I,R)
    I_normalised = I / normalisation
    coeff, var_matrix = curve_fit(gauss, R, I_normalised, guess)
    return coeff, var_matrix, normalisation

i = 3
centres = np.load("../data/detection_sweep/variable_aperture_r5/centres_" + str(i) + "sigma.npy")
galactic_intensities = np.load("../data/detection_sweep/variable_aperture_r5/galactic_intensities_" + str(i) + "sigma.npy")
background_intensities = np.load("../data/detection_sweep/variable_aperture_r5/background_intensities_" + str(i) + "sigma.npy")
img_data = np.load("../data/detection_sweep/variable_aperture_r5/final_img_" + str(i) + "sigma.npy")

hdulist = fits.open("../data/fits/mosaic.fits")
data = hdulist[0].data
data = data.astype(np.float64)

full_catalogue, stars, galaxies, spread_galaxies, indeterminate = [], [], [], [], []
sersic_cov_threshold = 0.0625
#moffat_cov_threshold = 0.05
moffat_beta_threshold_upper, moffat_beta_threshold_lower = 2, 1

for i in range(len(centres)):
    if i % 100 == 0:
        print ("iter: ", i)
    is_sersic, is_moffat = False, False
    x, y, r = centres[i]
    current_data = data - background_intensities[i]         # Account for Local background for current bright spot
    
    # Obtain Sersic Parameters
    total_I, total_pix =  s_p.circle_intensity(current_data, (x,y), r)
    if s_p.half_light_radius(current_data, (x,y), r, total_I) == None:
        continue
    Re, Re_corrected = s_p.half_light_radius(current_data, (x,y), r, total_I)               # Half light radius
    I_surface = [s_p.surface_intensity(current_data, (x, y), i) for i in range(r+1)]        # Surface intensities (ydata)
    Ie = I_surface[Re]                                                                      # Half light Intensity
    n_opt, sersic_cov = sersic_fit(Ie, Re, I_surface[2:], range(2, r + 1), 1)               # Obtain fitted sersic index n
    
    # Obtain Moffat Parameters
    beta_opt, moffat_cov = moffat_fit([i for i in range(r+1)], r, 0, I_surface[0],  I_surface, 1)
    moffat_cov = moffat_cov.flatten()
    current_moffat_cov_threshold = beta_opt[0] * 0.1
    
    if sersic_cov < sersic_cov_threshold:
        is_sersic = True
    
    if  beta_opt > 1 and moffat_cov < current_moffat_cov_threshold:
        is_moffat = True
        
    if not is_sersic and is_moffat:
        # Moffat Fitted => is a star
        full_catalogue.append([i, "2", beta_opt[0], moffat_cov[0]])
        stars.append([i, beta_opt[0]])
    else:
        # Galaxy
        full_catalogue.append([i, "1", n_opt[0], sersic_cov[0]])
        galaxies.append([i, n_opt[0], sersic_cov[0], beta_opt[0], moffat_cov[0]])

print ("The number of galaxies are: ", len(galaxies))
print ("The number of stars are: ", len(stars))

#spread_galaxies = np.array(spread_galaxies)
#spread_indices = spread_galaxies[:,0]
#spread_beta = spread_galaxies[:,3]
#n_indices = spread_galaxies[:,1]
#
#for i in range(len(spread_beta)):
#    current_beta = spread_beta[i]
#    multiplier = (1 + 0.06 * current_beta)
#    n_indices[i] *= multiplier
#
#bins = np.arange(0.5, 10.5, 0.5)
#plt.figure()
#plt.hist(n_indices, bins)
