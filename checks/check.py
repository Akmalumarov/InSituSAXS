import json
from pathlib import Path

folder = Path(__file__).resolve().parents[1] / "results/synthetic/Pq"

expected = {"spheres_1": 15.3, "long_cylinders_1": 12.0}
passed = True

for model, R_true in expected.items():
    try:
        with open(folder / f"{model}.json") as f:
            R_fit = json.load(f)["parameters"]["radius_nm"]

        error = abs(R_fit / R_true - 1)
        ok = error <= 0.10
        print(f"{model}: {'PASS' if ok else 'FAIL'} | error = {error:.1%}")
        passed &= ok
    except (OSError, KeyError, TypeError, ValueError) as e:
        print(f"{model}: FAIL ({e})")
        passed = False

raise SystemExit(0 if passed else 1)