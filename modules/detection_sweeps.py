#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  6 04:04:05 2019

@author: ricktjwong
"""

import numpy as np
from astropy.io import fits
import math
import time
import peak_detection as pk

hdulist = fits.open("../data/fits/masked.fits")
data = hdulist[0].data
data = data.astype(np.float64)

# 5 sigma - 3482
# 3 sigma - 3456
mu, sigma = 3418.5636925, 12.58387389

def sweep_fixed_aperture():
    """
    Conduct a sweep of different fixed aperture sizes from r = 5 to r = 10
    """
    for i in range(5, 11):
        start = time.time()
        threshold = math.ceil(mu + 3 * sigma)
        galaxy_count = pk.GalaxyCount(data, threshold, initial_r = i,
                                      aperture_type = "fixed")
        galaxy_count.count_galaxies()
        end = time.time()
        print("time: ")
        print(end - start)
        print(len(galaxy_count.background_intensities))
        print(len(galaxy_count.galactic_intensities))
        np.save("centres_"+str(i)+"radius", galaxy_count.centres)
        np.save("background_intensities_"+str(i)+"radius", galaxy_count.background_intensities)
        np.save("galactic_intensities_"+str(i)+"radius", galaxy_count.galactic_intensities)
        np.save("final_img_"+str(i)+"radius", galaxy_count.data)


def sweep_variable_aperture():
    """
    Conduct a sweep of different background sigma values for the variable
    aperture detection algorithm, with galaxy pixel accepted as 1sigma to
    10sigma away from the background mean
    """
    for i in range(3, 4):
        start = time.time()
        threshold = math.ceil(mu + i * sigma)
        galaxy_count = pk.GalaxyCount(data, threshold, initial_r=5)
        galaxy_count.count_galaxies()
        end = time.time()
        print("time: ")
        print(end - start)
        print(len(galaxy_count.background_intensities))
        print(len(galaxy_count.galactic_intensities))
        np.save("centres_"+str(i)+"sigma", galaxy_count.centres)
        np.save("background_intensities_"+str(i)+"sigma", galaxy_count.background_intensities)
        np.save("galactic_intensities_"+str(i)+"sigma", galaxy_count.galactic_intensities)
        np.save("final_img_"+str(i)+"sigma", galaxy_count.data)

sweep_variable_aperture()
