"""Временный анализ PDB: парсинг, граф связей, выделение частей молекулы."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np

PDB = Path("data/Zalygin_DESY/molecule.pdb")

ELEM = {
    "C": 6, "H": 1, "N": 7, "O": 8, "P": 15, "S": 16,
}


def parse_pdb(path: Path):
    atoms = []  # (idx, name, elem, x, y, z)
    bonds = []
    for line in path.read_text().splitlines():
        if line.startswith(("HETATM", "ATOM")):
            idx = int(line[6:11])
            name = line[12:16].strip()
            elem = line[76:78].strip() or name[0]
            x, y, z = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            atoms.append((idx, name, elem, x, y, z))
        elif line.startswith("CONECT"):
            parts = line.split()
            i = int(parts[1])
            for j in parts[2:]:
                bonds.append((i, int(j)))
    return atoms, bonds


def build_graph(atoms, bonds):
    adj = defaultdict(set)
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    return adj


def components(adj, atoms):
    idxs = [a[0] for a in atoms]
    seen = set()
    comps = []
    for start in idxs:
        if start in seen:
            continue
        stack = [start]
        comp = []
        while stack:
            n = stack.pop()
            if n in seen:
                continue
            seen.add(n)
            comp.append(n)
            stack.extend(adj[n] - seen)
        comps.append(comp)
    return comps


def formula(comp, atoms):
    counts = defaultdict(int)
    for i in comp:
        elem = next(a[2] for a in atoms if a[0] == i)
        counts[elem] += 1
    return dict(counts)


def electrons(counts):
    return sum(ELEM[e] * n for e, n in counts.items())


def main():
    atoms, bonds = parse_pdb(PDB)
    adj = build_graph(atoms, bonds)
    comps = components(adj, atoms)
    print(f"atoms={len(atoms)} bonds={len(bonds)} components={len(comps)}")
    for ci, comp in enumerate(sorted(comps, key=len, reverse=True)):
        f = formula(comp, atoms)
        ne = electrons(f)
        # bounding box
        coords = np.array([[a[3], a[4], a[5]] for a in atoms if a[0] in comp])
        lo, hi = coords.min(0), coords.max(0)
        span = hi - lo
        print(f"\n--- component {ci}: {len(comp)} atoms, {ne} e- ---")
        print(f"  formula: {f}")
        print(f"  bbox span (nm): {np.round(span, 3)}")
        print(f"  atoms: {sorted(comp)[:5]}...{sorted(comp)[-5:]}")


if __name__ == "__main__":
    main()