##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.0
## Maintainer: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
##################################################

from math import sqrt
from skimage.feature import blob_dog, blob_log, blob_doh
import matplotlib.pyplot as plt


def find_peaks(data, t):
    blobs_log = blob_log(data, min_sigma=4, max_sigma=30, num_sigma=27, threshold=t)
    # Compute radii in the 3rd column.
    blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)
    return blobs_log


def plot_blobs(data, blobs):
    fig, ax = plt.subplots()
    #coords_radius = []
    #
    for blob in blobs:
        ax.set_title("LOG")
        ax.imshow(data, interpolation='nearest')
        y, x, r = blob
    #    coords_radius.append([x, y, r])
        c = plt.Circle((x, y), r, color='red', linewidth=2, fill=False)
        ax.add_patch(c)
    
    plt.show()
