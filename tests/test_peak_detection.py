##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.1
## Maintainers: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
## Description: Unit tests of peak detection algo
## on single galaxies and galaxies of various
## shapes, as well as a small section of the data
##################################################

import sys
sys.path.append("../")
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import utils.plot as pt
import modules.peak_detection as pk
import time

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 16
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.figsize'] = 5, 5
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.default'] = 'regular'

hdulist = fits.open("../data/fits/masked.fits")
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
    data = data[901:1101, 901:1101]
#    data = data[3200:3230, 2445:2475]
#    data = data[101:301, 101:301]
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.imshow(data)
    galaxy_count = pk.GalaxyCount(data, 3482)
    galaxy_count.count_galaxies()
    assert(len(galaxy_count.background_intensities) == 10)
    assert(len(galaxy_count.galactic_intensities) == 10)
    print(len(galaxy_count.background_intensities))
    print(galaxy_count.galactic_intensities)
    print(len(galaxy_count.galactic_intensities))
    print(galaxy_count.centres)    
    centres = np.array(galaxy_count.centres)
    pt.mark_detected_objects(centres, ax)
    plt.xticks(np.arange(0, 200, 25))
    plt.yticks(np.arange(0, 200, 25))
    plt.savefig("unmasked.pdf", dpi=3000)    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.imshow(galaxy_count.data)
    pt.mark_detected_objects(centres, ax)
    plt.xticks(np.arange(0, 200, 25))
    plt.yticks(np.arange(0, 200, 25)) 
    plt.savefig("masked.pdf", dpi=3000)    

start = time.time()
test_small_section(data)
end = time.time()
print("time: ")
print(end - start)

# mu, sigma = 3418.5636925, 12.58387389
# 5 sigma - 3482
# 3 sigma - 3456

#mu, sigma = 3418.5636925, 12.58387389
#import math

#for i in range(1, 11):
#    start = time.time()
#    threshold = math.ceil(mu + i * sigma)
#    galaxy_count = pk.GalaxyCount(data,threshold)
#    galaxy_count.count_galaxies()
#    end = time.time()
#    print("time: ")
#    print(end - start)
#    print(len(galaxy_count.background_intensities))
#    print(len(galaxy_count.galactic_intensities))
#    np.save("centres_"+str(i)+"sigma", galaxy_count.centres)
#    np.save("background_intensities_"+str(i)+"sigma", galaxy_count.background_intensities)
#    np.save("galactic_intensities_"+str(i)+"sigma", galaxy_count.galactic_intensities)
#    np.save("final_img_"+str(i)+"sigma", galaxy_count.data)

#plt.figure()
#plt.imshow(data)
#start = time.time()
#galaxy_count = pk.GalaxyCount(data, 3482)
#galaxy_count.count_galaxies()
#end = time.time()
#print("time: ")
#print(end - start)
#print(len(galaxy_count.background_intensities))
#print(len(galaxy_count.galactic_intensities))
#plt.figure()
#plt.imshow(galaxy_count.data)

