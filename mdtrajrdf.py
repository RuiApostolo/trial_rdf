#  import os
import argparse
#  import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rdf

try:
    import matplotlib
    matplotlib.use('agg')  # no interactive plotting, only save figures
    import pylab
    have_matplotlib = True
except ImportError:
    have_matplotlib = False

parser = argparse.ArgumentParser(
    description='Generate the RDF g(f) of a lammpstrj file')
parser.add_argument(
    '--trajfile',
    '-t',
    help='Path to the lammpstrj file',
    required=True)
parser.add_argument(
    '--datafile',
    '-d',
    help='Path to the pdb file',
    required=True)
parser.add_argument(
    '--pairselect1',
    '-a1',
    help='Atom selection query for first part of pair',
    required=True)
parser.add_argument(
    '--pairselect2',
    '-a2',
    help='Atom selection query for first part of pair',
    required=True)
parser.add_argument(
    '--outdir',
    '-dir',
    help='Directory to output files into',
    default='.')
parser.add_argument(
    '--outname',
    '-o',
    help='Name for Output Files',
    default='rdf')

args = parser.parse_args()

u = mda.Universe(
    args.datafile,
    args.trajfile,
    atom_style="id resid type charge x y z",
    traj_atom_style="id type x y z vx vy vz fx fy fz",
    format='LAMMPSDUMP',
    topology_format='DATA',
    )

print("Trajectory and Topology Loaded")
print(u)

sel1 = u.select_atoms(args.pairselect1)
sel2 = u.select_atoms(args.pairselect2)

ave_rdf = rdf.InterRDF(
    sel1,
    sel2,
    nbins=100,
    range=(0.0, 10.0),
    )

ave_rdf.run()

print("RDF run")

outfile = args.outdir+'/'+args.outname+'.dat'
with open(outfile, 'w') as output:
    for ave_rdf.bins, ave_rdf.rdf in zip(ave_rdf.bins, ave_rdf.rdf):
        output.write(
            "{ave_rdf.bins:8.3f} \t {ave_rdf.rdf:8.3f}\n".format(
                **vars()))

print("g(r) data written to {outfile!r}".format(**vars()))

if have_matplotlib:
    matplotlib.rc('font', size=14)
    matplotlib.rc('figure', figsize=(5, 4))
    pylab.clf()
    pylab.plot(ave_rdf.bins, ave_rdf.rdf, linewidth=3)
    pylab.xlabel(r"distance $r$ in $\AA$")
    pylab.ylabel(r"radial distribution function $g(r)$")
    pylab.savefig(args.outdir+'/'+args.outname+'.pdf')
    pylab.savefig(args.outdir+'/'+args.outname+'.png')
