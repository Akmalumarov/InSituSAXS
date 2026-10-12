"""Запуск фита форм-факторов для синтетических SAXS-кривых.

Обрабатывает обе кривые из data/synthetic/Pq:
- spheres_1.txt      -> модель "sphere"
- long_cylinders_1.txt -> модель "cylinder" (длинный цилиндр, только радиус)

Сохраняет в results/synthetic/Pq:
- <name>.json — параметры фита (radius_nm, polydispersity, scale, background)
- <name>.png  — сравнительный график exp vs fit

Запуск: python3 src/Fitting/run_fit_synthetic.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from Fitting.fit_form_factors import fit_from_file, fit_result_to_dict

DATA_DIR = Path("data/synthetic/Pq")
OUT_DIR = Path("results/synthetic/Pq")

# модель для каждой кривой
MODELS = {
    "spheres_1.txt": "sphere",
    "long_cylinders_1.txt": "cylinder",
}

# ограничение области фита (q в нм^-1). Цилиндр: только осцилляции до фона.
Q_RANGE = {
    "spheres_1.txt": None,  # сфера фитируется по всему диапазону
    "long_cylinders_1.txt": (0.0, 0.11),  # осцилляции заканчиваются ~0.1 нм^-1
}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for fname, model in MODELS.items():
        data_path = DATA_DIR / fname
        name = data_path.stem

        result = fit_from_file(data_path, model, pd=0.1, q_range=Q_RANGE[fname])

        # JSON с параметрами
        out_json = OUT_DIR / f"{name}.json"
        out_json.write_text(
            json.dumps(fit_result_to_dict(result), indent=2, ensure_ascii=False)
        )

        # Сравнительный график
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.loglog(result.q_nm, result.I_exp, "o", ms=3, label="experiment")
        ax.loglog(result.q_nm, result.I_fit, "-", lw=1.5, label=f"fit ({model})")
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