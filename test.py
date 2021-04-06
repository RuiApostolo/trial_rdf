import tempfile, os, sys
import argparse

parser = argparse.ArgumentParser(description='Generate the RDF g(f) of a lammpstrj file')
parser.add_argument('--trajfile', '-t', help='Path to the lammpstrj file', required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--pdbfile', '-pdb', help='Path to the pdb file')
group.add_argument('--dumpfile', '-d', help='Path to the dump file')
parser.add_argument('--pairselect1', '-a1', help='Atom selection query for first part of pair', required=True)
parser.add_argument('--pairselect2', '-a2', help='Atom selection query for first part of pair', required=True)
parser.add_argument('--outdir', '-dir', help='Directory to output files into', default='./output')
parser.add_argument('--outname', '-file', help='Name for Output Files', default='rdf')

args = parser.parse_args()


print args.dumpfile
print args.pdbfile
