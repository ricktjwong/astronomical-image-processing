##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.0
## Maintainers: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
## Description: Profile fitting for galactic
## objects
##################################################

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
    b = 2 * n ** (-1/3)
    return Ie * np.exp(-b * ( (R/Re)**(1/n) - 1) )


def gauss(R, sigma):
    return 1 / (np.sqrt(np.pi * 2) * sigma) * \
            np.exp(-(R * R) / (2 * sigma * sigma))


def moffat(r, R, background, I0, beta):
    return background + I0 / (1 + (r * r) / (R * R)) ** beta


def moffat_fit(r, R, background, I0, I, guess):
    def moffat_beta(r, beta):
        return background + I0/(1 + (r*r)/(R*R) ) ** beta
    coeff, var_matrix = curve_fit(moffat_beta, r, I, guess)
    return coeff, var_matrix


def sersic_fit(Ie, Re, I_data, R_data, guess):
    R_data, I_data = np.array(R_data), np.array(I_data)
    def sersic_n(R_data, n):
         b = 2*n**(-1/3)
         return Ie * np.exp(-b * ((R_data / Re) ** (1/n) - 1))
    coeff, var_matrix = curve_fit(f = sersic_n, xdata = R_data,
                                  ydata = I_data, p0 = guess)
    return coeff, var_matrix


def gauss_fit(R,I,guess):
    normalisation = np.trapz(I,R)
    I_normalised = I / normalisation
    coeff, var_matrix = curve_fit(gauss, R, I_normalised, guess)
    return coeff, var_matrix, normalisation


centres = np.load("../data/detection/centres_5sigma.npy")
galactic_intensities = np.load("../data/detection/galactic_intensities_5sigma.npy")
background_intensities = np.load("../data/detection/background_intensities_5sigma.npy")
img_data = np.load('../data/detection/final_img_5sigma.npy')

hdulist = fits.open("../data/fits/masked.fits")
data = hdulist[0].data
data = data.astype(np.float64)
data = data[901:1101, 901:1101]
centres = np.array([[133, 26, 11], [65, 122, 12], [78, 180, 10], [116, 168, 6], [194, 103, 5], [52, 107, 5], [91, 180, 5], [80, 165, 5], [35, 6, 5], [33, 84, 5]])

full_catalogue, stars, galaxies, indeterminate = [], [], [], []
sersic_cov_threshold = 0.01
moffat_cov_threshold = 0.05
moffat_beta_threshold_upper, moffat_beta_threshold_lower = 2, 1

for i in range(len(centres)):
    if i % 50 == 0:
        print ("iter: ", i)
    is_sersic, is_moffat = False, False
    x, y, r = centres[i]
    current_data = data - background_intensities[i]         # Account for Local background for current bright spot
    
    # Obtain Sersic Parameters
    total_I, total_pix = s_p.circle_intensity(current_data, (x, y), r)
    Re, Re_corrected = s_p.half_light_radius(current_data, (x, y), r, total_I)           # Half light radius
    I_surface = [s_p.surface_intensity(current_data, (x, y), i) for i in range(r+1)]    # Surface intensities (ydata)
    Ie = I_surface[Re]                                                                  # Half light Intensity
    n_opt, sersic_cov = sersic_fit(Ie, Re, I_surface[2:], range(2,r+1), 1)              # Obtain fitted sersic index n
    
    # Obtain Moffat Parameters
    beta_opt, moffat_cov = moffat_fit([i for i in range(r+1)], r, 0, I_surface[0],  I_surface, 1)
    
    if sersic_cov < sersic_cov_threshold:
        is_sersic = True
        
    if beta_opt < moffat_beta_threshold_upper and beta_opt > moffat_beta_threshold_lower:
        is_moffat = True
        
    if not is_sersic and not is_moffat:     # Indeterminate
        full_catalogue.append([i,"0",n_opt,sersic_cov,beta_opt,moffat_cov])
        indeterminate.append([i,n_opt,sersic_cov,beta_opt,moffat_cov])        
        
    if is_sersic and is_moffat:             # Galaxy
        full_catalogue.append([i,"1",n_opt,sersic_cov,beta_opt,moffat_cov])
        galaxies.append([i,n_opt,sersic_cov,beta_opt,moffat_cov])
    
    if is_sersic and not is_moffat:         # Sersic Fitted => is a galaxy
        full_catalogue.append([i,"1",n_opt,sersic_cov])
        galaxies.append([i,n_opt])
        
    if not is_sersic and is_moffat:         # Moffat Fitted => is a star
        full_catalogue.append([i,"2",beta_opt,moffat_cov])
        stars.append([i,beta_opt])

print ("The number of indeterminates are: ", len(indeterminate))
print ("The number of galaxies are: ", len(galaxies))
print ("The number of stars are: ", len(stars))    

plt.figure(1)
plt.imshow(data)
plt.show()    
