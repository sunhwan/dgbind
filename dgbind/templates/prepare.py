import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import argparse
import numpy as np
import random

def weighted_choice(weights):
    totals = np.cumsum(weights)
    rnd = random.random() * totals[-1]
    for i, total in enumerate(totals):
        if rnd < total:
            return i

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dcd', metavar='dcdfile', help='input DCD filename')
    parser.add_argument('pdb', metavar='pdbfile', help='output PDB filename')
    parser.add_argument('target', metavar='target', help='target RMSD value')
    args = parser.parse_args()

    universe = mda.Universe("{{ psffile }}", args.dcdfile)
    ref = mda.Universe("{{ psffile }}", "{{ pdbfile }}")
    R = RMSD(universe, ref, select="{{ selection }}")
    R.run()

    kbt = 1.0e-3 * 8.314472 / 4.1868 * {{ temperature }}
    target = float(args.target)
    P = np.exp(-1/kbt*0.5*{{ force }}*(R.rmsd[:,2] - target)**2)
    P = P / np.sum(P)
    idx = weighted_choice(P)
    ts = universe.trajectory[idx]
    mda.Writer(args.pdbfile, bonds=False).write(universe)