#!/usr/bin/python
import pylab as plt
import numpy as np

F_fast = [[],[],[]]
F_therm = [[],[],[]]
F_tot = [[],[],[]]

fl = open("meshtal")
parse = False
for line in fl:
    if(parse):
        line = line.strip()
        line = line.split()
        if(line[0] == "6.250E-07"):
            F_therm[0].append(float(line[1][:5]))
            F_therm[1].append(float(line[2].replace('-','E-')))
            F_therm[2].append(float(line[3]))

        elif(line[0] == "1.000E+02"):
            F_fast[0].append(float(line[1][:5]))
            F_fast[1].append(float(line[2].replace('-','E-')))
            F_fast[2].append(float(line[3]))

        elif(line[0] == "Total"):
            F_tot[0].append(float(line[1][:5]))
            F_tot[1].append(float(line[2].replace('-','E-')))
            F_tot[2].append(float(line[3]))

    elif (("Energy" in line) and ("Rel Error" in line)):
        parse = True
fl.close()

F_max = 0.0
for i in range(len(F_fast[1])):
    if(F_fast[1][i] > F_max):
        F_max = F_fast[1][i]
    #if(F_therm[1][i] > F_max):
    #    F_max = F_therm[1][i]

for i in range(len(F_fast[1])):
    F_fast[1][i] = F_fast[1][i]/F_max
    #F_therm[1][i] = F_therm[1][i]/F_max

F_fast = np.array(F_fast)
#F_therm = np.array(F_therm)
F_tot = np.array(F_tot)

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size

font_size = plt.rcParams["font.size"]
font_size = 15
plt.rcParams["font.size"] = font_size

plt.errorbar(F_fast[0], F_fast[1], yerr=F_fast[1]*F_fast[2], c='b', label="Flux Rapide")
#plt.errorbar(F_therm[0], F_therm[1], yerr=F_therm[1]*F_therm[2], c='orange', label="Flux Thermique")
#plt.errorbar(F_tot[0], F_tot[1], yerr=F_tot[1]*F_tot[2], c='r', label="Total")
plt.xlabel("Position [cm]")
plt.ylabel("Flux")
#plt.legend()
plt.tight_layout()
plt.show()

plt.plot(F_fast[0], F_fast[1]*F_fast[2], c='b', label="Erreur Flux Rapide")
#plt.errorbar(F_therm[0], F_therm[1]*F_therm[2], c='orange', label="Erreur Flux Thermique")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.legend()
plt.xlabel("Position [cm]")
plt.ylabel("Flux Erreur")
plt.tight_layout()
plt.show()
