from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

hdulist = fits.open("../data/mosaic.fits")
data = hdulist[0].data
print("data:")
print(data)
print (type(data))
print (len(data))
print (data[0])

def get_brightest(nparray):
    max = 0
    for i in nparray:
        if i > max:
            max = i
    return max

x = np.array([[1,1,1,1,1], [1,1,1,1,2]])

for i in range(200):
    print(get_brightest(data[i]))