#!/usr/bin/python3

import matplotlib.pyplot as plt

L = 2.0
NFOMBIN = 100
dx = L / NFOMBIN
x = []
for i in range(100):
    x.append(i*dx + 0.5*dx)

XS = {}
fl = open("Coll_Densities.txt")
line_count = 0;
for line in fl:
    line = line.strip()
    line = line.split(",")
    if(len(line) < 2):
        continue

    if(line[0] == "#XS"):
        xsname = line[1]
        XS[xsname] = TM = {} # Dic for data of all tracking methods
    elif(line[0] == "#TM"):
        tmname = line[1]
        XS[xsname][tmname] = [[], [], []]
        line_count = 0
    else:
        for elem in line:
            XS[xsname][tmname][line_count].append(float(elem))
        line_count += 1
fl.close()

#print(XS.keys())
#print(XS["BroadGaussian"].keys())
#print(XS["BroadGaussian"]["DT"])

xsnames = ["LI", "LD", "EI", "ED", "SG", "BG"]
tmnames = ["DS", "DT", "MDT", "NWDT", "MNWDT", "BT", "MBT"]
#tmnames = ["DT", "MDT", "NWDT", "MNWDT"]

for xs in xsnames:
    # Plot coll. density
    for tm in tmnames:
        plt.plot(x,XS[xs][tm][0],label=tm)
    plt.title(xs+" Collision Density")
    plt.legend()
    plt.show()

    # Plot FOM
    for tm in tmnames:
        plt.plot(x,XS[xs][tm][2],label=tm)
    plt.title(xs+" FOM")
    plt.legend()
    plt.show()
