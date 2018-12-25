#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 14:27:33 2018

@author: ricktjwong
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

z = 25.3

hdulist = fits.open("../data/masked.fits")
data = hdulist[0].data
data = data.astype(np.float64)
    
data = data[101:801, 101:801]

plt.figure(1)
plt.imshow(data)

class GalaxyCount:
    def __init__(self, data):
        self.data = data
        self.rows = data.shape[0]
        self.cols = data.shape[1]
        # Set threshold as 5 standard deviations away from mean background
        self.threshold = 3485
        # Stores the total intensity of galaxies of radius >= 5 pixels
        self.galactic_intensities = []
        # Stores the mean intensity of background corresponding to the region
        # where the galaxy was detected
        self.background_intensities = []
        
    def get_brightest(self):
        """
        Get the brightest pixel in the image data and set the indices
        """
        self.max_flattened = np.argmax(self.data)
        self.max_idx_r = int(self.max_flattened / self.cols)
        self.max_idx_c = self.max_flattened % self.cols    
        
    def set_aperture(self):
        """
        Defines an aperture within which we find pixels corresponding to
        galaxies or background
        """
        self.box_x1 = self.max_idx_c - self.r
        self.box_x2 = self.max_idx_c + self.r
        self.box_y1 = self.max_idx_r - self.r    
        self.box_y2 = self.max_idx_r + self.r
        
    def build_catalogue(self):
        """
        Loop through the dimensions of the aperture and append to the local
        background or galaxy list
        """
        for i in range(self.box_y1, self.box_y2 + 1):
            for j in range(self.box_x1, self.box_x2 + 1):
                if self.data[i][j] < self.threshold:
                    self.background_intensity.append(data[i][j])
                else:
                    self.galaxy_intensity.append(data[i][j])
                    self.data[i][j] = 0                    

    def get_catalogue(self):
        """
        Does radius and invalid index checks before building the catalogue
        """
        # Set lists to store background and galaxy intensities. These are for
        # individual galaxies within an aperture
        self.set_aperture()
        # Check if the aperture is out of range - if so, mask the point
        if self.box_x1 < 1 or self.box_x2 > self.cols - 1 or self.box_y1 < 1 \
        or self.box_y2 > self.rows - 1:
            self.data[self.max_idx_r][self.max_idx_c] = 0
            self.get_brightest()
            self.n += 1
        # If the aperture is within range of the image data, build catalogue
        else:
            self.background_intensity = []
            self.galaxy_intensity = []
            self.build_catalogue()
            # If the length of background intensity is < galaxy pixels,
            # increase the radius by 1 pixel and continue appending
            while len(self.background_intensity) < len(self.galaxy_intensity):
                # Reset the individual galaxy and corresponding background list
                self.background_intensity = []
                self.galaxy_intensity = []
                self.r += 1
                # Set the aperture based on new radius
                self.set_aperture()
                # Build catalog based on new aperture
                self.build_catalogue()
            # Only accept a galaxy if it spans at least 15 pixels in area
            if len(self.galaxy_intensity) > 15:
                self.galactic_intensities += [sum(self.galaxy_intensity)]
                self.background_intensities += [np.mean(self.background_intensity)]
            self.n += 1

    def count_galaxies(self):
        """
        Repeat till all galaxies are found
        """
        self.n = 0
        self.galactic_intensities = []
        self.background_intensities = []
        #while (max_flattened != 0):
        while (self.n < 200):
            # Reset aperture radius to 5
            self.r = 5
            # Get new brightest soot in the image
            self.get_brightest()
            # Do checks and build catalogue
            self.get_catalogue()
    
galaxy_count = GalaxyCount(data)
galaxy_count.count_galaxies()

print(len(galaxy_count.background_intensities))
print(len(galaxy_count.galactic_intensities))

plt.figure(2)
plt.imshow(galaxy_count.data)

def save_to_fits(data):
    hdu = fits.PrimaryHDU(data)
    hdul = fits.HDUList([hdu])
    hdul.writeto('complete.fits')
