#!/usr/bin/python3

import matplotlib.pyplot as plt

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

xsname = "LinearlyDecreasing"

#plt.plot(XS[xsname]["DS"][0])
plt.plot(XS[xsname]["DT"][2])
#plt.plot(XS[xsname]["MDT"][0])
plt.plot(XS[xsname]["NWDT"][2])
plt.plot(XS[xsname]["MNWDT"][2])
plt.plot(XS[xsname]["BT"][2])
#plt.plot(XS[xsname]["MBT"][0])
plt.show()
