##################################################
## Author: Rick Wong and Daniel Seah
## Version: 1.0.1
## Maintainers: ricktjwong and danielsrq
## Email: rtw16@ic.ac.uk and drs16@ic.ac.uk
## Description: General shared plotting module
##################################################

import matplotlib.pyplot as plt


def mark_detected_objects(centres, ax):
    for i in range(len(centres)):
        circle = plt.Circle((centres[i][0], centres[i][1]), centres[i][2],
                            color='r',  fill=False)
        ax.add_artist(circle)
    plt.scatter(centres[:,0], centres[:,1], c='r', s=5)
