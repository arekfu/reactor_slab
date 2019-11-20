#!/usr/bin/python3

import matplotlib.pyplot as plt

flux = []
nvc = []
fl = open("data.csv")

line = fl.readline()
line = line.split(",")
for el in line:
    flux.append(float(el))

line = fl.readline()
line = line.split(",")
for el in line:
    nvc.append(float(el))

plt.plot(flux)
plt.show()

plt.plot(nvc)
plt.show()
