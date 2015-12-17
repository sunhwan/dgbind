import MDAnalysis as mda
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

    universe = mda.Universe("{{ psffile }}", args.dcd)
    ref = mda.Universe("{{ psffile }}", "{{ pdbfile }}")
{% if jobtype == 'RMSD' %}
    from MDAnalysis.analysis.rms import RMSD
    R = RMSD(universe, ref, select="{{ selection }}")
    R.run()
    obs = R.rmsd[:,2]
{% elif jobtype == 'Angles' %}
    from MDAnalysis.core.Timeseries import TimeseriesCollection, CenterOfGeometry
    import numpy.linalg as la
    collection = TimeseriesCollection()
    refatoms = []
    {% for key in angle -%}
    refatoms.append(universe.select_atoms("{{ refatoms[key] }}"))
    {% endfor -%}
    for refatom in refatoms:
        collection.addTimeseries(CenterOfGeometry(refatom))
    collection.compute(universe.trajectory)

    r = []
    natoms = len(collection)
    for i in range(natoms-1):
        r.append((np.array(collection[i+1]) - np.array(collection[i])).T)
    if natoms == 3:
        v = r[0][:,0]*r[1][:,0] + r[0][:,1]*r[1][:,1] + r[0][:,2]*r[1][:,2]
        theta = 180 - np.rad2deg(np.arccos(v / la.norm(r[0], axis=1) / la.norm(r[1], axis=1)))
    if natoms == 4:
        v1 = np.cross(r[0], r[1])
        v2 = np.cross(r[1], r[2])
        n1 = v1 / la.norm(v1, axis=1)[:,np.newaxis]
        n2 = v2 / la.norm(v2, axis=1)[:,np.newaxis]
        x = n1[:,0]*n2[:,0] + n1[:,1]*n2[:,1] + n1[:,2]*n2[:,2]
        
        m1 = np.cross(n1, r[1]/la.norm(r[1], axis=1)[:,np.newaxis])
        y = m1[:,0]*n2[:,0] + m1[:,1]*n2[:,1] + m1[:,2]*n2[:,2]
        theta = -np.rad2deg(np.arctan2(y, x))
    obs = theta
{% elif jobtype == 'Distance' %}
    from MDAnalysis.core.Timeseries import TimeseriesCollection, CenterOfGeometry
    collection = TimeseriesCollection()
    p1 = universe.select_atoms("{{ refatoms.p1 }}")
    l1 = universe.select_atoms("{{ refatoms.l1 }}")
    collection.addTimeseries(CenterOfGeometry(p1))
    collection.addTimeseries(CenterOfGeometry(l1))
    collection.compute(universe.trajectory)
    dist = np.sqrt(np.sum((np.array(collection[0]) - np.array(collection[1]))**2, axis=0))
    obs = dist
{% endif %}

    kbt = 1.0e-3 * 8.314472 / 4.1868 * {{ temperature }}
    target = float(args.target)
    P = np.exp(-1/kbt*0.5*{{ force }}*(obs - target)**2)
    P = P / np.sum(P)
    idx = weighted_choice(P)
    ts = universe.trajectory[idx]
    mda.Writer(args.pdb, bonds=False).write(universe)