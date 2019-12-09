#!/usr/bin/python3

# Et = x + y + z + 1
# Es = 0.7*Et(r) , Ec = 0.3*Et(r)

import numpy as np
import openmc

nxbins = 30
nybins = 30 
nzbins = 30

L = 2.0

dx = L/nxbins
dy = L/nybins
dz = L/nzbins

groups = openmc.mgxs.EnergyGroups(np.logspace(-5,7,3))

# Make all cross sections
XSDATAS = []
for i in range(nxbins):
    for j in range(nybins):
        for k in range(nzbins):
            xsname = str(i)+"."+str(j)+"."+str(k)
            xs = openmc.XSdata(xsname, groups)
            xs.order = 0
            
            # Get position at center
            x = i*dx + 0.5*dx
            y = j*dy + 0.5*dy
            z = k*dz + 0.5*dz
            # Get XS values
            Et = x + y + z + 1.0
            Ec = 0.3*Et
            Es = 0.7*Et
            # assign XSs
            xs.set_total([Et,Et], temperature=294.)
            xs.set_absorption([Ec,Ec], temperature=294.)
            xs.set_nu_fission([0.0,0.0], temperature=294.)
            sct_matrx = [[[0.5*Es, 0.5*Es],
                          [0.0, Es]]]
            sct_matrx = np.array(sct_matrx)
            sct_matrx = np.rollaxis(sct_matrx, 0, 3)
            xs.set_scatter_matrix(sct_matrx, temperature=294.)
            # Add xs to XSDATAS
            XSDATAS.append(xs)

mg_xs_file = openmc.MGXSLibrary(groups)
for i in range(len(XSDATAS)):
    mg_xs_file.add_xsdata(XSDATAS[i])
mg_xs_file.export_to_hdf5('mgxsd.h5')

# Make Materials
materials = {}
for i in range(nxbins):
    for j in range(nybins):
        for k in range(nzbins):
            xsname = str(i)+"."+str(j)+"."+str(k)
            materials[xsname] = openmc.Material(name=xsname)
            materials[xsname].set_density('macro', 1.)
            materials[xsname].add_macroscopic(xsname)

mat_file = openmc.Materials(materials.values())
mat_file.cross_sections = 'mgxsd.h5'
mat_file.export_to_xml()

# Make Geometry
nxsurfs = nxbins + 1
nysurfs = nybins + 1
nzsurfs = nzbins + 1

XPLANES = []
YPLANES = []
ZPLANES = []
for i in range(nxsurfs):
    if (i == 0) or (i == nxsurfs - 1):
        XPLANES.append(openmc.XPlane(x0 = i*dx, boundary_type='vacuum'))
    else:
        XPLANES.append(openmc.XPlane(x0 = i*dx))

for j in range(nysurfs):
    if (j == 0) or (j == nysurfs - 1):
        YPLANES.append(openmc.YPlane(y0 = j*dy, boundary_type='vacuum'))
    else:
        YPLANES.append(openmc.YPlane(y0 = j*dy))

for k in range(nzsurfs):
    if (k == 0) or (k == nzsurfs - 1):
        ZPLANES.append(openmc.ZPlane(z0 = k*dz, boundary_type='vacuum'))
    else:
        ZPLANES.append(openmc.ZPlane(z0 = k*dz))

CELLS = []
for i in range(nxbins):
    for j in range(nybins):
        for k in range(nzbins):
            xsname = str(i)+"."+str(j)+"."+str(k)
            CELLS.append(openmc.Cell())
            CELLS[-1].fill = materials[xsname]
            CELLS[-1].region = +XPLANES[i] & -XPLANES[i+1] & +YPLANES[j] & -YPLANES[j+1] & +ZPLANES[k] & -ZPLANES[k+1]

universe = openmc.Universe()
universe.add_cells(CELLS)

geom_file = openmc.Geometry(universe)
geom_file.export_to_xml()

# Settings
particles = 1000000
batches = 50
inactive = 0

sett_file = openmc.Settings()
sett_file.batches = batches
sett_file.inactive = inactive
sett_file.particles = particles
sett_file.verbosity = 4
sett_file.energy_mode = 'multi-group'
sett_file.run_mode = 'fixed source'
point = openmc.stats.Point(xyz=(1.,1.,1.))
sett_file.source = openmc.Source(space=point)
sett_file.export_to_xml()








