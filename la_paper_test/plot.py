#!/usr/bin/python3

import matplotlib.pyplot as plt

# List of bins tested
bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# 0 = col_std, 1 = esc_std, 2 = xs_eval, 3 = FOM_coll, 4 = FOM_esc
XS_types = ['LI', 'LD', 'EI', 'ED', 'SG', 'BG']
Trak_types = ['DS', 'DT', 'MDT', 'NWDT', 'MNWDT']

data = {}
# Initialize empty data
for xs in XS_types:
    for trak in Trak_types:
        data_key = xs+"-"+trak
        data[data_key] = [[], [], [], [], []]

# Get data from files
for bin in bins:
    fname = str(bin) + "-bins.txt"
    fl = open(fname)

    # read to lines to get rid of header
    line = fl.readline()
    line = fl.readline()

    xs = -1
    trk = -1
    cnt = 0
    parse = True
    while parse:
        line = fl.readline()
        if not line: break
        line = line.strip()
        line = line.replace(",","")
        if("------" in line): 
            xs += 1
            line = fl.readline()
        elif("Direct" in line): trk = 0
        elif(line == "Delta Tracking"): trk = 1
        elif("Meshed Delta" in line): trk = 2
        elif(line == "Negative Weight Delta Tracking"): trk = 3
        elif("Meshed Neg" in line): trk = 4
        elif(line == ""): continue
        else:
            key = XS_types[xs]+"-"+Trak_types[trk]
            line = line.split(" ")
            data[key][0].append(float(line[4]))
            data[key][1].append(float(line[9]))
            line = fl.readline()
            line = line.strip()
            line = line.replace(",","")
            line = line.split()
            data[key][2].append(float(line[-1]))
            line = fl.readline()
            line = line.strip()
            line = line.replace(",","")
            line = line.split()
            data[key][3].append(float(line[2]))
            data[key][4].append(float(line[5]))

    fl.close()

# Make plots

# MNWDT nbins
for xs in XS_types:
    key = xs+"-MNWDT"
    plt.plot(bins, data[key][2], label=xs)
plt.legend()
plt.xlabel("Bins")
plt.ylabel("Average Number of XS Evaluations")
plt.title("Meshed Negative Weighted Delta Tracking")
plt.show()

for xs in XS_types:
    key = xs+"-MNWDT"
    plt.plot(bins, data[key][0], label=xs)
plt.legend()
plt.xlabel("Bins")
plt.ylabel("Collision STD")
plt.title("Meshed Negative Weighted Delta Tracking")
plt.show()

for xs in XS_types:
    key = xs+"-MNWDT"
    plt.plot(bins, data[key][1], label=xs)
plt.legend()
plt.xlabel("Bins")
plt.ylabel("Transmision STD")
plt.title("Meshed Negative Weighted Delta Tracking")
plt.show()

for xs in XS_types:
    key = xs+"-MNWDT"
    plt.plot(bins, data[key][3], label=xs)
plt.legend()
plt.xlabel("Bins")
plt.ylabel("Collision FOM")
plt.title("Meshed Negative Weighted Delta Tracking")
plt.show()

for xs in XS_types:
    key = xs+"-MNWDT"
    plt.plot(bins, data[key][4], label=xs)
plt.legend()
plt.xlabel("Bins")
plt.ylabel("Transmision FOM")
plt.title("Meshed Negative Weighted Delta Tracking")
plt.show()

# MDT nbins
for xs in XS_types:
    key = xs+"-MDT"
    plt.plot(bins, data[key][0], label=xs)
plt.legend()
plt.xlabel("Bins")
plt.ylabel("Collision STD")
plt.title("Meshed Delta Tracking")
plt.show()

for xs in XS_types:
    key = xs+"-MDT"
    plt.plot(bins, data[key][3], label=xs)
plt.legend()
plt.xlabel("Bins")
plt.ylabel("Collision FOM")
plt.title("Meshed Delta Tracking")
plt.show()

