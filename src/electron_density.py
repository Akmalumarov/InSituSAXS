"""Электронная плотность слоёв молекулы OQ8 (система Залыгина).

Молекула — двойной хвост (две C18-цепи) + полярная голова (фосфат, N/O/S).
Для core-shell цилиндра нужны средние электронные плотности:
  - ядро (core)  = гидрофобные хвосты;
  - корона (corona) = полярная голова (+ растворитель в реальной системе).

Метод объёмов:
  - хвост: групповые объёмы CH2=27, CH3=54 A^3 (стандарт для алкильных цепей);
  - голова: парциальный молярный объём = V_mol - V_tail, где V_mol оценён
    из ван-дер-ваальсова объёма с коэффициентом упаковки f_pack.

rho = N_e / V  (e/A^3).
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np

PDB = Path("data/Zalygin_DESY/molecule.pdb")

ELEM = {"C": 6, "H": 1, "N": 7, "O": 8, "P": 15, "S": 16}
# ван-дер-ваальсовы радиусы, A
VDW_R = {"C": 1.7, "H": 1.2, "N": 1.55, "O": 1.52, "P": 1.8, "S": 1.8}
# групповые объёмы алкильной цепи, A^3
V_CH2, V_CH3 = 27.0, 54.0
# коэффициент упаковки vdW -> парциальный молярный объём (для головы)
F_PACK = 1.72

# две C18-цепи (индексы атомов углерода в PDB)
TAIL_C = set(range(2, 20)) | set(range(57, 75))


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


def build_graph(bonds):
    adj = defaultdict(set)
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    return adj


def split_tail_head(atoms, adj):
    """Хвост = две C18-цепи + их H; голова = всё остальное."""
    tail = set(TAIL_C)
    for a in atoms:
        if a[2] == "H" and any(m in TAIL_C for m in adj[a[0]]):
            tail.add(a[0])
    head = {a[0] for a in atoms} - tail
    return tail, head


def composition(atoms, idxs):
    counts = defaultdict(int)
    for a in atoms:
        if a[0] in idxs:
            counts[a[2]] += 1
    return dict(counts)


def electrons(comp):
    return sum(ELEM[e] * n for e, n in comp.items())


def vdw_volume(atoms, idxs, grid=160):
    """Ван-дер-ваальсов объём набора атомов (метод сетки)."""
    sel = [a for a in atoms if a[0] in idxs]
    coords = np.array([[a[3], a[4], a[5]] for a in sel])
    radii = np.array([VDW_R[a[2]] for a in sel])
    lo = coords.min(0) - radii.max()
    hi = coords.max(0) + radii.max()
    g = [np.linspace(lo[d], hi[d], grid) for d in range(3)]
    X, Y, Z = np.meshgrid(*g, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    inside = np.zeros(len(pts), bool)
    for c, r in zip(coords, radii):
        inside |= np.linalg.norm(pts - c, axis=1) <= r
    return inside.sum() * ((hi - lo).prod() / grid**3)


def tail_group_volume(comp):
    """Объём хвоста из групповых объёмов: 2 цепи, каждая C18H33 = 17 CH2 + 1 CH3."""
    n_ch2 = 2 * 17
    n_ch3 = 2 * 1
    return n_ch2 * V_CH2 + n_ch3 * V_CH3


def main():
    atoms, bonds = parse_pdb(PDB)
    adj = build_graph(bonds)
    tail, head = split_tail_head(atoms, adj)

    comp_tail = composition(atoms, tail)
    comp_head = composition(atoms, head)
    ne_tail = electrons(comp_tail)
    ne_head = electrons(comp_head)

    # объёмы
    v_tail = tail_group_volume(comp_tail)
    v_mol_vdw = vdw_volume(atoms, set(adj) | {a[0] for a in atoms})
    v_mol = v_mol_vdw * F_PACK          # парциальный молярный объём молекулы
    v_head = v_mol - v_tail             # парциальный объём головы

    rho_tail = ne_tail / v_tail
    rho_head = ne_head / v_head
    rho_mol = (ne_tail + ne_head) / v_mol

    print(f"молекула OQ8: {len(atoms)} атомов, {ne_tail + ne_head} e-")
    print(f"  хвост (core):  {comp_tail}  N_e={ne_tail}")
    print(f"  голова (corona): {comp_head}  N_e={ne_head}")
    print(f"  V_tail (групповые) = {v_tail:.0f} A^3")
    print(f"  V_mol (vdW={v_mol_vdw:.0f} * {F_PACK}) = {v_mol:.0f} A^3")
    print(f"  V_head = {v_head:.0f} A^3")
    print(f"\n  rho_core  = {rho_tail:.4f} e/A^3")
    print(f"  rho_corona(голова) = {rho_head:.4f} e/A^3")
    print(f"  rho_molecule = {rho_mol:.4f} e/A^3")
    print(f"\n  справочно: вода rho=0.334 e/A^3; углеводород ~0.28 e/A^3")
    print(f"  корона в реальной системе включает растворитель: rho_corona = "
          f"(N_e,head + N_e,solv) / V_corona")


if __name__ == "__main__":
    main()