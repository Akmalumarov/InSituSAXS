"""Запуск фита форм-факторов для синтетических SAXS-кривых.

Обрабатывает обе кривые из data/synthetic/Pq:
- spheres_1.txt      -> модель "sphere"
- long_cylinders_1.txt -> модель "cylinder" (длинный цилиндр, только радиус)

Область фита для цилиндра определяется по данным:
- нижняя граница = 0.5 * q_first_minima (первый минимум форм-фактора);
- верхняя граница = q, где кривая выходит на фоновое плато
  (интенсивность приближается к минимальному значению хвоста).

Сохраняет в results/synthetic/Pq:
- <name>.json — параметры фита (radius_nm, polydispersity, scale, background)
- <name>.png  — сравнительный график exp vs fit с границами области фита

Запуск: python3 src/Fitting/run_fit_synthetic.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from Fitting.fit_form_factors import fit_form_factor, fit_result_to_dict

DATA_DIR = Path("data/synthetic/Pq")
OUT_DIR = Path("results/synthetic/Pq")

# модель для каждой кривой
MODELS = {
    "spheres_1.txt": "sphere",
    "long_cylinders_1.txt": "cylinder",
}


def find_first_minimum(q_nm: np.ndarray, I: np.ndarray) -> float:
    """Первый локальный минимум лог-интенсивности (переход убывание->возрастание)."""
    logI = np.log(I + 1e-12)
    d = np.diff(logI)
    sign = np.sign(d)
    idx = np.where((sign[:-1] < 0) & (sign[1:] > 0))[0]
    if len(idx) > 0:
        return float(q_nm[idx[0] + 1])
    return float(q_nm[-1])


def find_plateau_start(q_nm: np.ndarray, I: np.ndarray) -> float:
    """q, где кривая выходит на фоновое плато.

    Плато — хвост кривой, где интенсивность приближается к фоновому значению
    (минимуму кривой). Считаем, что плато начинается там, где I впервые
    опускается ниже 2× минимального значения (фон ~ I_min).
    """
    I_min = float(np.min(I))
    thr = 2.0 * I_min
    below = np.where(I < thr)[0]
    if len(below) == 0:
        return float(q_nm[-1])
    return float(q_nm[below[0]])


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for fname, model in MODELS.items():
        data_path = DATA_DIR / fname
        name = data_path.stem

        data = np.loadtxt(data_path, comments="#", skiprows=1)
        q_ang = data[:, 0]
        I_exp = data[:, 1]
        q_nm = q_ang * 10.0  # Å^-1 -> нм^-1

        # Определяем область фита по данным
        q_range = None
        if model == "cylinder":
            q_first_min = find_first_minimum(q_nm, I_exp)
            q_lo = 0.5 * q_first_min
            q_hi = find_plateau_start(q_nm, I_exp)
            q_range = (q_lo, q_hi)
            print(f"[{name}] первый минимум q={q_first_min:.4f} нм^-1")
            print(f"[{name}] область фита: q in [{q_lo:.4f}, {q_hi:.4f}] нм^-1")

        result = fit_form_factor(q_ang, I_exp, model, pd=0.1, q_range=q_range)

        # JSON с параметрами
        out_json = OUT_DIR / f"{name}.json"
        out_json.write_text(
            json.dumps(fit_result_to_dict(result), indent=2, ensure_ascii=False)
        )

        # Сравнительный график с границами области фита
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.loglog(result.q_nm, result.I_exp, "o", ms=3, label="experiment")
        ax.loglog(result.q_nm, result.I_fit, "-", lw=1.5, label=f"fit ({model})")
        if q_range is not None:
            qlo, qhi = q_range
            ax.axvline(qlo, color="tab:red", ls="--", lw=1, label=f"q_min=0.5·q1st ({qlo:.3f})")
            ax.axvline(qhi, color="tab:green", ls="--", lw=1, label=f"q_plateau ({qhi:.3f})")
        ax.set_xlabel("q (nm$^{-1}$)")
        ax.set_ylabel("I(q)")
        ax.set_title(f"{name} — {model}, R = {result.radius_nm:.3f} nm")
        ax.legend()
        fig.tight_layout()
        fig.savefig(OUT_DIR / f"{name}.png", dpi=150)
        plt.close(fig)

        print(f"[{name}] model={model} success={result.success}")
        print(f"  radius_nm={result.radius_nm:.4f}  pd={result.polydispersity}")
        print(f"  scale={result.scale:.4e}  background={result.background:.4e}")
        print(f"  -> {out_json}")


if __name__ == "__main__":
    main()