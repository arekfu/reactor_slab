#!/usr/bin/python3

import matplotlib.pyplot as plt

print(" 1) DT   2) MDT   3) NWDT   4)MNWDT   5) CT   6) MCT   7)IMCT")
track = int(input(" Select Tracking Method => "))
if track == 1:
    tm = "dt"
elif track == 2:
    tm = "mdt"
elif track == 3:
    tm = "nwdt"
elif track == 4:
    tm = "mnwdt"
elif track == 5:
    tm = "ct"
elif track == 6:
    tm = "mct"
elif track == 7:
    tm = "imct"
else:
    exit()

L = 2.0
NFOMBIN = 100
dx = L / NFOMBIN
x = []
for i in range(NFOMBIN):
    x.append(i*dx + 0.5*dx)

XS = {}
fl = open("li_"+tm+"_coll_density.txt")
#fl = open("CONTROL.txt")
line_count = 0; 
# 0 = real avg, 1 = real std, 2 = real FOM, 3 = all avg, 4 = all std, 5 = all FOM
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
        XS[xsname][tmname] = [[], [], [], [], [], []]
        line_count = 0
    else:
        for elem in line:
            XS[xsname][tmname][line_count].append(float(elem))
        line_count += 1
fl.close()

xsnames = ["LI"]
tmnames = [tm.upper()]

for xs in xsnames:
    # Plot coll. density
    for tm in tmnames:
        plt.plot(x,XS[xs][tm][0], label=tm)
        plt.plot(x,XS[xs][tm][3], label=tm+" all")
        #plt.errorbar(x,XS[xs][tm][0], yerr=XS[xs][tm][1],label=tm)
    plt.title(xs+" Collision Density")
    plt.legend()
    plt.show()

    # Plot FOM
    for tm in tmnames:
        plt.plot(x,XS[xs][tm][2],label=tm)
        plt.plot(x,XS[xs][tm][5],label=tm+" all")
    plt.title(xs+" FOM")
    plt.legend()
    plt.show()
