#!/usr/bin/python3

import matplotlib.pyplot as plt

FOM_DT_CNT = []
FOM_DT_T = []
FOM_NDT_CNT = []
FOM_NDT_T = []

x = []

fl = open("FOM_output.txt")

cnt = -1
for line in fl:
    line = line.strip()
    line = line.split(",")
    if cnt == -1:
        p = float(line[0])
        q = float(line[1])
    for fom in line:
        if cnt == 0:
            FOM_DT_CNT.append(float(fom))
        elif cnt == 1:
            FOM_DT_T.append(float(fom))
        elif cnt == 2:
            FOM_NDT_CNT.append(float(fom))
        elif cnt == 3:
            FOM_NDT_T.append(float(fom))
    cnt += 1

L = 2.0
dx = 2.0 / len(FOM_DT_T)
for i in range(len(FOM_DT_T)):
    x.append(i*dx + 0.5*dx)

plt.plot(x, FOM_DT_CNT, label="Delta Tracking")
plt.plot(x, FOM_NDT_CNT, label="Negative Delta Tracking")
plt.title("Average XS Evals FOM - p = "+str(p)+", q = "+str(q))
plt.legend()
plt.show()

plt.plot(x, FOM_DT_T, label="Delta Tracking")
plt.plot(x, FOM_NDT_T, label="Negative Delta Tracking")
plt.title("Time FOM - p = "+str(p)+", q = "+str(q))
plt.legend()
plt.show()
