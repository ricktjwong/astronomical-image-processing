##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.0
## Maintainer: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
##################################################

import numpy as np


def circle_intensity(array, centre, radius):
    if type(radius) == float:
        raise TypeError("Radius has to be int")
    if radius < 0:
        raise ValueError("Radius has to be at least 1")
    max_y, min_y = centre[1] + radius, centre[1] - radius
    max_x, min_x = centre[0] + radius, centre[0] - radius
    cumulative_sum = 0
    total = 0
    for i in range(min_y, max_y + 1):
        for j in range(min_x, max_x + 1):
            if ((i - centre[1]) **2  + (j - centre[0])**2 ) <= (radius * radius):
                total += 1
#                print("[" + str(i) + "," + str(j) + "]")
                cumulative_sum += array[i][j]                
    return cumulative_sum, total


def half_light_radius(array, centre, full_radius, total_intensity):
    half_intensity = total_intensity * 0.5
    previous_intensity = 0
    for i in range(1, full_radius):
        current_intensity = circle_intensity(array, centre, i)[0]
        if current_intensity >= half_intensity:
            correction = (current_intensity - half_intensity) / (current_intensity - previous_intensity)    # Linear correction
            return (i, i - correction)          # returns integer radius and linearly corrected radius 
        previous_intensity = current_intensity


def surface_intensity(data, centre, radius):
    if radius == 0:
        return data[centre[1], centre[0]]
    intensity_1, pixels_1 = circle_intensity(data, centre, radius)
    intensity_2, pixels_2 = circle_intensity(data, centre, radius - 1)
    return (intensity_1 - intensity_2) / (pixels_1 - pixels_2)


def I(R, I_e, R_e, n):
    b = 2 * (n ** (-1/3))
    return I_e * np.exp(-b * (((R/R_e) ** (1/n)) - 1))
