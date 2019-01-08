#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 10:25:17 2019

@author: ricktjwong
"""

import numpy as np
import matplotlib.pyplot as plt

galaxies = np.load("../data/classification/galaxies-0.0625.npy")
stars = np.load("../data/classification/stars-0.0625.npy")

indices = galaxies[:,0]
beta = galaxies[:,3]
n_indices = galaxies[:,1]

plt.figure()
bins = np.arange(0.5, 10.5, 0.5)
plt.hist(n_indices, bins)

for i in range(len(beta)):
    current_beta = beta[i]
    multiplier = (1 + 0.06 * current_beta)
    n_indices[i] *= multiplier

bins = np.arange(0.5, 10.5, 0.5)
plt.figure()
plt.hist(n_indices, bins, edgecolor='black', color='b')
