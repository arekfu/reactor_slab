#!/usr/bin/python3

import pylab as plt

fl = open("terminal_output")
fl.readline()
fl.readline()
fl.readline()
fl.readline()

gen = []
keff = []
kavg = [[], []]
H = []
FOM = []

ignore = True
for line in fl:
    if(("---" in line) or ("Gen" in line)):
        continue

    line = line.split()
    if(len(line) > 3):
        ignore = False

    elif(len(line) == 0):
        break

    gen.append(float(line[0]))
    keff.append(float(line[1]))
    H.append(float(line[2]))
    if(ignore == False):
        kavg[0].append(float(line[4]))
        kavg[1].append(float(line[5]))
        FOM.append(float(line[6]))
fl.close()

ngens = gen[-1]
nignored = ngens - len(kavg[0])

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size

font_size = plt.rcParams["font.size"]
font_size = 15
plt.rcParams["font.size"] = font_size

plt.plot(gen,keff, label="Keff")
plt.xlabel("Génération")
plt.tight_layout()
plt.show()

#plt.errorbar(gen[int(nignored):], kavg[0], yerr=kavg[1], c='r')
plt.plot(gen[int(nignored):], kavg[0], c='r')
plt.xlabel("Génération")
plt.tight_layout()
plt.show()

plt.plot(gen, H, c="g")
plt.show()

plt.plot(gen[int(nignored):], FOM, c="orange")
plt.show()
