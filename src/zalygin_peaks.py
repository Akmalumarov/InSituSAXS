"""Анализ SAXS-кривых образцов Залыгина (DESY): поиск пиков, определение фазы, графики."""
from __future__ import annotations

import numpy as np
from pathlib import Path
from scipy.signal import savgol_filter, find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = Path("data/Zalygin_DESY")
OUT = Path("results/Zalygin/Desy/графики")

# Соотношения q/q1 для лиотропных фаз
PHASES = {
    "hexagonal HII": [1, np.sqrt(3), 2, np.sqrt(7), 3],
    "lamellar":      [1, 2, 3, 4, 5],
    "cubic Im3m":    [np.sqrt(2), np.sqrt(4), np.sqrt(6), np.sqrt(8), np.sqrt(10)],
    "cubic Pn3m":    [np.sqrt(2), np.sqrt(3), np.sqrt(4), np.sqrt(6), np.sqrt(8)],
    "cubic Ia3d":    [np.sqrt(6), np.sqrt(8), np.sqrt(14), np.sqrt(16), np.sqrt(20)],
    "FCC":           [np.sqrt(3), np.sqrt(4), np.sqrt(8), np.sqrt(11), np.sqrt(12)],
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
    """Сглаживание + поиск пиков в диапазоне q."""
    mask = (q >= qmin) & (q <= qmax)
    qs, Is = q[mask], I[mask]
    if len(qs) < window:
        return []
    Ism = savgol_filter(Is, window, 2)
    prom = prom_frac * Ism.max()
    peaks, _ = find_peaks(Ism, prominence=prom)
    return [(qs[p], Ism[p]) for p in peaks]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, name, color in zip(axes, ("red", "green"), ("crimson", "seagreen")):
        q, I, err = load_curve(name)
        ax.loglog(q, I, color=color, lw=1.2, label=name)
        # пики в области структурного фактора
        peaks = find_peaks_smooth(q, I, 0.02, 0.6, window=15, prom_frac=0.03)
        for qp, Ip in peaks:
            ax.axvline(qp, color="k", ls="--", lw=0.8, alpha=0.6)
            ax.annotate(f"{qp:.3f}", (qp, Ip), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=8)
        ax.set_xlabel("q (Å⁻¹)")
        ax.set_ylabel("I (a.u.)")
        ax.set_title(f"{name} — {len(peaks)} пиков")
        ax.legend()
        print(f"=== {name} === пики (q, I):")
        for qp, Ip in peaks:
            print(f"  q={qp:.5f}  I={Ip:.3e}")
    fig.tight_layout()
    fig.savefig(OUT / "zalygin_curves_loglog.png", dpi=150)
    print("saved:", OUT / "zalygin_curves_loglog.png")


if __name__ == "__main__":
    main()