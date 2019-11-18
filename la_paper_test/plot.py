#!/usr/bin/python3

import matplotlib.pyplot as plt

nbins = [1,2,3,4,5]

LI = [2.53, 2.28, 2.64, 3.13, 3.67]
LD = [1.42, 1.48, 1.68, 1.94, 2.21]
EI = [7.10, 3.58, 3.36, 3.68, 4.16]
ED = [2.24, 2.28, 2.82, 3.47, 4.17]
SG = [2.47, 3.44, 3.38, 4.60, 5.10]
BG = [1.25, 1.80, 2.32, 2.87, 3.43]

plt.plot(nbins, LI, label="Linearly Increasing")
plt.plot(nbins, LD, label="Linearly Decreasing")
plt.plot(nbins, EI, label="Exponentially Increasing")
plt.plot(nbins, ED, label="Exponentially Decreasing")
plt.plot(nbins, SG, label="Sharp Gaussian")
plt.plot(nbins, BG, label="Broad Gaussian")
plt.legend()
plt.xlabel("Number of Delta Tracking Bins")
plt.ylabel("Average Evaluations per Collision")
plt.show()
