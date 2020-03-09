#!/usr/bin/env python
import tempfile, os, sys
import argparse
# import numpy as np
from dump import dump
from pdbfile import pdbfile
# import mdtraj as md
import mdanalysis as mda
path = os.environ["LAMMPS_PYTHON_TOOLS"]
sys.path.append(path)

try:
    import matplotlib
    matplotlib.use('agg')  # no interactive plotting, only save figures
    import pylab
    have_matplotlib = True
except ImportError:
    have_matplotlib = False

# receive parameters from command line
parser = argparse.ArgumentParser(description='Generate the RDF g(f) of a lammpstrj file')
parser.add_argument('--trajfile', '-t', help='Path to the lammpstrj file', required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--pdbfile', '-pdb', help='Path to the pdb file')
group.add_argument('--dumpfile', '-d', help='Path to the dump file')
parser.add_argument('--pairselect1', '-a1', help='Atom selection query for first part of pair', required=True)
parser.add_argument('--pairselect2', '-a2', help='Atom selection query for first part of pair', required=True)
parser.add_argument('--outdir', '-dir', help='Directory to output files into', default='./output')
parser.add_argument('--outname', '-f', help='Name for Output Files', default='rdf')
parser.add_argument('--dumporder', '-o', help='Columns #s for ID,type,x,y,z (usually 1,2,3,4,5)', default='1,2,3,4,5')

# Parse arguments
args = parser.parse_args()

# Split dumporder into list
dol = args.dumporder.split(',')

# use pizza.py to automatically convert to pdb
# requires LAMMPS python tool subfolder to be set as an env variable:
# export LAMMPS_PYTHON_TOOLS="~/LAMMPS/20170331/tools/python/pizza"

# use pizza.py's dump routine to get pdb file
# d = dump(args.dumpfile, 0)
# d.map(int(dol[0]), "id", int(dol[1]), "type", int(dol[2]), "x", int(dol[3]), "y", int(dol[4]), "z")
# time = d.next()
# d.aselect.all()

# p = pdbfile(d)
# p.one(args.outname)

# modified:
# try:
#     dumpfile
# except NameError:
#     traj = md.load(args.trajfile, top=args.pdbfile)
#     print("Trajectory and PDB Topology Loaded")
# else:
#     dump2pdb(args.dumpfile, 1, 2, 3, 4, 5, file.txt)

# original
traj = md.load(args.trajfile, top=args.pdbfile)
print("Trajectory and Topology Loaded")

# print("done")
pairs = traj.top.select_pairs(args.pairselect1, args.pairselect2)
radii, rdf = md.geometry.rdf.compute_rdf(traj, pairs, r_range=(0.0, 10.0), bin_width=0.1)
print("Computed")

outfile = args.outdir+'/'+args.outname+'.dat'
with open(outfile, 'w') as output:
    for radius, gofr in zip(radii, rdf):
        output.write("{radius:8.3f} \t {gofr:8.3f}\n".format(**vars()))
print ("g(r) data written to {outfile!r}".format(**vars()))

if have_matplotlib:
    matplotlib.rc('font', size=14)
    matplotlib.rc('figure', figsize=(5, 4))
    pylab.clf()
    pylab.plot(radii, rdf, linewidth=3)
    pylab.xlabel(r"distance $r$ in $\AA$")
    pylab.ylabel(r"radial distribution function $g(r)$")
    pylab.savefig(args.outdir+'/'+args.outname+'.pdf')
    pylab.savefig(args.outdir+'/'+args.outname+'.png')
