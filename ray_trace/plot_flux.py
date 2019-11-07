#!/usr/bin/python3
import matplotlib
from matplotlib import pyplot as plt

#matplotlib.rcParams['figure.dpi'] = 200

x = []
F = [[],[]]
F_err = [[],[]]

# Read flux data file
fl = open("ray_trace_flux.csv")
for i in range(5):
    line = fl.readline()
    line = line.split(",")
    if (i == 0):
        # Add values to x array
        for j in range(len(line)):
            x.append(float(line[j]))
    elif (i == 1):
        # Fast flux
        for j in range(len(line)):
            F[0].append(float(line[j]))
    elif (i == 2):
        # Thermal flux
        for j in range(len(line)):
            F[1].append(float(line[j]))
    elif (i == 3):
        # Fast flux error
        for j in range(len(line)):
            F_err[0].append(float(line[j]))
    elif (i == 4):
        # Thermal flux error
        for j in range(len(line)):
            F_err[1].append(float(line[j]))

fl.close()

F_max = 0.0
# Normalize flux and error
for i in range(len(F[0])):
    if(F[0][i] > F_max):
        F_max = F[0][i]
    if(F[1][i] > F_max):
        F_max = F[1][i]
for i in range(len(F[0])):
    F[0][i] = F[0][i]/F_max
    F[1][i] = F[1][i]/F_max
    F_err[0][i] = F_err[0][i]/F_max
    F_err[1][i] = F_err[1][i]/F_max

# Set font sizes and such. Dont remember what
# the first bit does quite frankly mdr
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size

font_size = plt.rcParams["font.size"]
font_size = 15
plt.rcParams["font.size"] = font_size

# Error bar cap size in points
cpsize = 1

# Plot flux
plt.errorbar(x,F[0], yerr=F_err[0], c='blue', label="Flux Rapide", capsize=cpsize)
plt.errorbar(x,F[1], yerr=F_err[1], c='orange', label="Flux Thermique", capsize=cpsize)
plt.xlabel("Position [cm]")
plt.ylabel("Flux")
plt.tight_layout()
plt.show()

# plot flux error
plt.plot(x, F_err[0], c='blue', label="Erreur Flux Rapide")
plt.plot(x, F_err[1], c='orange', label="Erreur Flux Thermique")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel("Position [cm]")
plt.ylabel("Flux Erreur")
plt.tight_layout()
plt.show()
