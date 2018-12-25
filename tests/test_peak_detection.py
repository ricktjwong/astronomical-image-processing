#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 25 16:29:17 2018

@author: ricktjwong
"""

import sys
sys.path.append("../")
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import modules.peak_detection as pk

hdulist = fits.open("../data/masked.fits")
data = hdulist[0].data
data = data.astype(np.float64)


def generate_2D_gaussian(sigma):
    l, r = -10, 10
    x, y = np.meshgrid(np.linspace(l, r, 50), np.linspace(l, r, 50))
    d = np.sqrt(x * x + y * y)
    mu = 0.0
    # 2D gaussian distribution for a spherically symmetric source
    g = np.exp(-((d-mu)**2 / (2.0 * sigma**2)))
    return g


def test_single_symmetric_galaxy():
    """
    Generate a mock galaxy using a 2D gaussian distribution and assert that
    only one galaxy is detected
    """
    data = generate_2D_gaussian(2.0)
    plt.figure()
    plt.imshow(data)
    galaxy_count = pk.GalaxyCount(data, 0.01)
    galaxy_count.count_galaxies()
    assert(len(galaxy_count.background_intensities) == 1)
    assert(len(galaxy_count.galactic_intensities) == 1)
    plt.figure()
    plt.imshow(galaxy_count.data)


def test_single_symmetric_small_galaxy():
    """
    Generate a mock galaxy using a 2D gaussian distribution and assert that
    only one galaxy is detected
    """
    data = generate_2D_gaussian(0.2)
    plt.figure()
    plt.imshow(data)
    galaxy_count = pk.GalaxyCount(data, 0.01)
    galaxy_count.count_galaxies()
    print(galaxy_count.background_intensities)
    print(galaxy_count.centres)
    assert(len(galaxy_count.background_intensities) == 0)
    assert(len(galaxy_count.galactic_intensities) == 0)
    plt.figure()
    plt.imshow(galaxy_count.data)


def test_small_section(data):
    data = data[101:301, 101:301]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.imshow(data)
    galaxy_count = pk.GalaxyCount(data, 3485)
    galaxy_count.count_galaxies()
    assert(len(galaxy_count.background_intensities) == 10)
    assert(len(galaxy_count.galactic_intensities) == 10)
    print(len(galaxy_count.background_intensities))
    print(len(galaxy_count.galactic_intensities))
    print(galaxy_count.centres)    
    centres = np.array(galaxy_count.centres)
    plt.scatter(centres[:,0], centres[:,1], c='r', s=5)
    for i in range(len(centres)):
        circle = plt.Circle((centres[i][0], centres[i][1]), centres[i][2],
                            color='r',  fill=False)
        ax.add_artist(circle)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.imshow(galaxy_count.data)
    for i in range(len(centres)):
        circle = plt.Circle((centres[i][0], centres[i][1]), centres[i][2],
                            color='r',  fill=False)
        ax.add_artist(circle)
    plt.scatter(centres[:,0], centres[:,1], c='r', s=5)

test_small_section(data)
