#!/usr/bin/python3

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})

L = 2.0
NFOMBIN = 100
dx = L / NFOMBIN
x = []
for i in range(NFOMBIN):
    x.append(i*dx + 0.5*dx)

XS = {}
fl = open("Coll_Densities.txt")
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
        XS[xsname][tmname] = [[], [], [], []]
        line_count = 0
    else:
        for elem in line:
            XS[xsname][tmname][line_count].append(float(elem))
        line_count += 1
fl.close()

xsnames = ["LI","LD","EI","ED","SG","BG"]
#xsnames = ["LI"]
#tmnames = ["DT", "MDT", "NWDT", "MNWDT", "BT", "MBT", "PBT"]
tmnames = ["DT", "MDT", "NWDT", "MNWDT", "BT", "MBT", "IBT"]
tmnames = ["DT", "NWDT", "CT"]

x_li = [1.7,1.7]
y_li = [-1.0, 1.0]
x_bg_1 = [0.826864,0.826864]
x_bg_2 = [1.63314, 1.63314]
y_bg = [-1.0, 1.0]

for xs in xsnames:
    # Plot coll. density
    for tm in tmnames:
        plt.plot(x,XS[xs][tm][0], label=tm)
#        plt.plot(x,XS[xs][tm][2], label=tm+" all")
        #plt.errorbar(x,XS[xs][tm][0], yerr=XS[xs][tm][1],label=tm)
    plt.tight_layout()
    plt.title("Collision Density")
    plt.legend()
    plt.xlabel("Position")
    plt.show()

    # Plot FOM
    for tm in tmnames:
        plt.plot(x,XS[xs][tm][1],label=tm)
#        plt.plot(x,XS[xs][tm][3],label=tm+" all")
    if(xs == "LI"):
        plt.plot(x_li,y_li,c="black", lw=0.5)
        plt.ylim([-0.0007,0.0075])
        #plt.plot(x_bg_1,y_li,c="black", lw=0.5)
        #plt.plot(x_bg_2,y_li,c="black", lw=0.5)
    plt.tight_layout()
    plt.title("Figure of Merit")
    plt.legend()
    plt.xlabel("Position")
    plt.show()
