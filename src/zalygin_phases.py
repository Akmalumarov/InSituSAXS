"""Анализ SAXS-кривых образцов Залыгина (DESY). q в nm^-1.
Определение лиотропной фазы по соотношениям пиков + графики."""
from __future__ import annotations

import numpy as np
from pathlib import Path
from scipy.signal import savgol_filter, find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = Path("data/Zalygin_DESY")
OUT = Path("results/Zalygin/Desy/графики")

# Соотношения q/q1 для лиотропных фаз (первые 4-5 рефлексов)
PHASES = {
    "hexagonal HII": ([1, np.sqrt(3), 2, np.sqrt(7), 3], ["10", "11", "20", "21", "30"]),
    "lamellar":      ([1, 2, 3, 4, 5], ["001", "002", "003", "004", "005"]),
    "cubic Im3m":    ([np.sqrt(2), np.sqrt(4), np.sqrt(6), np.sqrt(8), np.sqrt(10)],
                      ["110", "200", "211", "220", "310"]),
    "cubic Pn3m":    ([np.sqrt(2), np.sqrt(3), np.sqrt(4), np.sqrt(6), np.sqrt(8)],
                      ["110", "111", "200", "211", "220"]),
    "cubic Ia3d":    ([np.sqrt(6), np.sqrt(8), np.sqrt(14), np.sqrt(16), np.sqrt(20)],
                      ["211", "220", "321", "400", "420"]),
    "FCC":           ([np.sqrt(3), np.sqrt(4), np.sqrt(8), np.sqrt(11), np.sqrt(12)],
                      ["111", "200", "220", "311", "222"]),
}


def load_curve(name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = []
    with open(DATA / f"{name}.dat") as f:
        for line in f:
            p = line.split()
            if len(p) == 3:
                rows.append([float(x) for x in p])
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2]


def find_peaks_smooth(q: np.ndarray, I: np.ndarray, qmin: float, qmax: float,
                      window: int, prom_frac: float) -> list[tuple[float, float]]:
    mask = (q >= qmin) & (q <= qmax)
    qs, Is = q[mask], I[mask]
    if len(qs) < window:
        return []
    Ism = savgol_filter(Is, window, 2)
    prom = prom_frac * Ism.max()
    peaks, _ = find_peaks(Ism, prominence=prom)
    return [(qs[p], Ism[p]) for p in peaks]


def match_phase(peaks: list[float], tol: float = 0.05) -> tuple[str, list[str], float]:
    """Определяет фазу по отношениям q/q1. Штраф за пропущенные рефлексы.
    Возвращает (имя, индексы, q1).

    Данные: пики из find_peaks_smooth (кривые data/Zalygin_DESY/red.dat, green.dat)
    Проверка: задача 1 (определение фазы) — checks/check.py, погрешность q ≤ 5%
    Используется в: src/zalygin_phases.py::main
    """
    q1 = peaks[0]
    ratios = [p / q1 for p in peaks]
    best, best_score, best_idx = None, 1e9, None
    for name, (ref, idx) in PHASES.items():
        refn = [r / ref[0] for r in ref]  # нормировка на первый рефлекс
        fit = sum(min(abs(r - r0) / r0 for r0 in refn) for r in ratios)
        # штраф за рефлексы фазы в диапазоне [1, max_ratio], не покрытые пиками
        maxr = max(ratios)
        missing = sum(1 for r0 in refn if 1 < r0 < maxr
                      and not any(abs(r - r0) / r0 < tol for r in ratios))
        score = fit + 0.1 * missing
        if score < best_score:
            best, best_score, best_idx = name, score, idx
    return best, best_idx, q1


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    summary = {}
    for ax, name, color in zip(axes, ("red", "green"), ("crimson", "seagreen")):
        q, I, err = load_curve(name)
        ax.loglog(q, I, color=color, lw=1.3, label=f"{name} (эксперимент)")
        peaks = find_peaks_smooth(q, I, 0.5, 1.5, window=15, prom_frac=0.05)
        qp = [p[0] for p in peaks]
        phase, idx, q1 = match_phase(qp)
        # индексы по ближайшему рефлексу (пики могут пропускать рефлексы)
        refn = [r / PHASES[phase][0][0] for r in PHASES[phase][0]]
        assigned = []
        for p in qp:
            r = p / q1
            j = min(range(len(refn)), key=lambda k: abs(r - refn[k]))
            assigned.append(PHASES[phase][1][j])
        # параметр решётки
        if phase == "hexagonal HII":
            a = 4 * np.pi / (np.sqrt(3) * q1)
            lat = f"a = {a:.2f} nm"
        elif phase == "FCC":
            a = 2 * np.pi * np.sqrt(3) / q1
            lat = f"a = {a:.2f} nm"
        else:
            lat = ""
        for (qpi, Ipi), mi in zip(peaks, assigned):
            ax.axvline(qpi, color="k", ls="--", lw=0.9, alpha=0.7)
            ax.annotate(f"({mi})\n{qpi:.3f}", (qpi, Ipi),
                        textcoords="offset points", xytext=(0, 10),
                        ha="center", fontsize=8, color="black")
        ax.set_xlabel("q (nm⁻¹)")
        ax.set_ylabel("I (a.u.)")
        ax.set_title(f"{name}: {phase}\n{lat}", fontsize=11)
        ax.legend(loc="upper right")
        summary[name] = {"phase": phase, "peaks": qp, "indices": idx, "lattice": lat}
        print(f"=== {name} === {phase} {lat}")
        for (qpi, Ipi), mi in zip(peaks, assigned):
            print(f"  ({mi}) q={qpi:.4f} nm^-1  I={Ipi:.3e}")
    fig.tight_layout()
    fig.savefig(OUT / "zalygin_phases.png", dpi=150)
    print("saved:", OUT / "zalygin_phases.png")
    np.save(OUT / "zalygin_peaks.npy", np.array(summary, dtype=object), allow_pickle=True)


if __name__ == "__main__":
    main()