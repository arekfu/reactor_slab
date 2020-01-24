#!/usr/bin/python3

import matplotlib.pyplot as plt

TMs = ["dt", "mdt", "nwdt", "mnwdt", "ct", "mct", "imct"]
x = [1,2,3,4,5,6,7]
lbls = ["DT", "MDT", "NWDT", "MNWDT", "CT", "MCT", "IMCT"]

k  = []
k_std = []
FOM = []

for tm in TMs:
    file = open("li_"+tm+".out")
    for line in file:
        if "Avg keff" in line:
            line = line.split()
            k.append(float(line[3]))
            k_std.append(float(line[5]))
            FOM.append(float(line[-1]))
    file.close()

plt.errorbar(x,k,yerr=k_std,linestyle="None", capsize=4, capthick=2, marker=".", ms=12)
plt.xticks(x,lbls)
plt.xlabel("Traking Method")
plt.ylabel("Keff")
plt.tight_layout()
plt.show()

plt.scatter(x, FOM)
plt.xticks(x,lbls)
plt.xlabel("Traking Method")
plt.ylabel("FOM")
plt.tight_layout()
plt.show()

