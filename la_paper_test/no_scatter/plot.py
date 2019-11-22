#!/usr/bin/python3

import matplotlib.pyplot as plt

nbins = [1,2,3,4,5,6,7,8,9,10]

DT = [2.53, 1.42, 7.10, 2.24, 2.49, 1.33]

LI = [2.39, 1.54, 1.30, 1.18, 1.12, 1.07, 1.04, 1.02, 1.00, 0.99]
LD = [1.28, 1.12, 1.05, 1.01, 0.98, 0.96, 0.95, 0.94, 0.93, 0.93]
EI = [7.03, 2.78, 1.90, 1.58, 1.42, 1.32, 1.25, 1.20, 1.17, 1.14]
ED = [1.52, 0.84, 0.63, 0.53, 0.48, 0.44, 0.42, 0.40, 0.39, 0.37]
SG = [1.55, 0.76, 0.53, 0.38, 0.52, 0.26, 0.23, 0.36, 0.17, 0.26]
BG = [1.02, 0.99, 0.89, 0.84, 0.81, 0.79, 0.78, 0.77, 0.76, 0.75]

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
