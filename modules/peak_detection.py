##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.1
## Maintainers: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
## Description: Peak detection algorithm which
## makes use of a variable 2D square aperture
##################################################

import numpy as np

z = 25.3

class GalaxyCount:
    def __init__(self, data, threshold, initial_r=5, aperture_type="variable"):
        self.data = data
        self.rows = data.shape[0]
        self.cols = data.shape[1]        
        # Set threshold as n standard deviations away from mean background
        self.threshold = threshold
        self.initial_r = initial_r
        self.aperture_type = aperture_type
        # Only accept a galaxy if it spans at least 15 pixels in area
        self.area_acceptance = 15
        # Stores the total intensity of galaxies of radius >= 5 pixels
        self.galactic_intensities = []
        # Stores the mean intensity of background corresponding to the region
        # where the galaxy was detected
        self.background_intensities = []
        self.centres = []
        
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
        update_data = self.data.copy()
        box_aperture = update_data[self.box_y1:self.box_y2 + 1,
                                   self.box_x1:self.box_x2 + 1]
        self.background_intensity = box_aperture[box_aperture <
                                                 self.threshold].copy()
        self.galaxy_intensity = box_aperture[box_aperture >=
                                             self.threshold].copy()
        # Mask the pixels if they are accepted as a galactic point
        box_aperture[box_aperture >= self.threshold] = 0
        return update_data

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
        # If the aperture is within range of the image data, build catalogue
        else:
            update_data = self.build_catalogue().copy()
            # If the length of background intensity is < galaxy pixels,
            # increase the radius by 1 pixel and continue appending
            if self.aperture_type == "variable":
                while len(self.background_intensity) < len(self.galaxy_intensity) \
                and len(self.galaxy_intensity) - len(self.background_intensity) != 1:
                    self.r += 1
                    # Set the aperture based on new radius
                    self.set_aperture()
                    # Build catalog based on new aperture
                    update_data = self.build_catalogue().copy()
            self.data = update_data.copy()
            if len(self.galaxy_intensity) >= self.area_acceptance:
                if self.aperture_type == "variable":
                    background_mean = np.mean(self.background_intensity)
                else:
                    background_mean = 3418.5636925
                self.galactic_intensities += \
                [sum(np.array(self.galaxy_intensity) - background_mean)]
                self.background_intensities += [background_mean]
                self.centres.append([self.max_idx_c, self.max_idx_r, self.r])

    def count_galaxies(self):
        """
        Repeat till all galaxies are found
        """
        self.n = 0
        self.galactic_intensities = []
        self.background_intensities = []
        #while (max_flattened != 0):
        while(True):
            # Reset aperture radius to initial start radius
            self.r = self.initial_r
            # Get new brightest soot in the image
            self.get_brightest()
            if self.data[self.max_idx_r][self.max_idx_c] < self.threshold:
                break
            if self.n % 1000 == 0:
                print(self.n)
                print(self.data[self.max_idx_r][self.max_idx_c])
            # Do checks and build catalogue
            self.get_catalogue()
            self.n += 1
