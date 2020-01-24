import matplotlib.pyplot as plt

file = open('Fission_Output.txt')

source = []
for line in file:
    source.append(float(line))
file.close()

max = 2.0
dx = 0.1
bins = []
nbins = int(max / dx)
for i in range(nbins):
    bins.append(0.0)

for p in source:
    bin = int(p/dx)
    if bin < nbins:
        bins[bin] += 1.0

plt.plot(bins)
plt.show()
