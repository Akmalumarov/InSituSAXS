"""Форм-факторы рассеяния для SAXS-моделей.

Соглашения:
- q в nm^-1, размеры в nm, электронные плотности в e-/nm^3.
- Контраст rho задаётся относительно растворителя (rho_solv=0 по умолчанию).
- Амплитуды F(q) возвращаются в единицах объёма частицы (V * контраст).
- Интенсивность P(q) = |F(q)|^2 / V (нормировка на объём, как в SAS-моделях).

Для цилиндров с произвольным радиальным профилем электронной плотности
используется численное интегрирование по радиусу (см. cylinder_radial_*).
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.integrate import quad
from scipy.special import j0, j1


# ---------------------------------------------------------------------------
# Вспомогательные устойчивые функции
# ---------------------------------------------------------------------------

def _sinc(x: np.ndarray) -> np.ndarray:
    """sin(x)/x, устойчиво при x -> 0.

    Данные: x (массив аргументов)
    Проверка: нет независимой проверки (вспомогательная)
    Используется в: src/FormFactors/form_factors.py::cylinder_amplitude, ::core_shell_cylinder_amplitude, ::cylinder_radial_amplitude
    """
    return np.sinc(np.asarray(x, dtype=float) / np.pi)


def _j1_over_x(x: np.ndarray) -> np.ndarray:
    """J1(x)/x, устойчиво при x -> 0 (предел = 1/2).

    Данные: x (массив аргументов)
    Проверка: нет независимой проверки (вспомогательная)
    Используется в: src/FormFactors/form_factors.py::infinite_cylinder_amplitude, ::cylinder_amplitude, ::core_shell_cylinder_amplitude, ::infinite_cylinder_intensity_pd, ::infinite_cylinder_intensity_pd_norm
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    out[small] = 0.5 - x[small] ** 2 / 16.0  # разложение J1(x)/x
    out[~small] = j1(x[~small]) / x[~small]
    return out


def _sphere_f(x: np.ndarray) -> np.ndarray:
    """f(x) = 3 (sin x - x cos x) / x^3, f(0) = 1.

    Данные: x (массив аргументов)
    Проверка: нет независимой проверки (вспомогательная)
    Используется в: src/FormFactors/form_factors.py::sphere_amplitude, ::core_shell_sphere_amplitude, ::sphere_intensity_pd, ::sphere_intensity_pd_norm
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    out[small] = 1.0 - x[small] ** 2 / 10.0  # разложение
    out[~small] = 3.0 * (np.sin(x[~small]) - x[~small] * np.cos(x[~small])) / x[~small] ** 3
    return out


def _orient_avg(amp_fn: Callable[[float], float], q: np.ndarray, V: float) -> np.ndarray:
    """Усреднение по ориентации: P(q) = (1/V) ∫_0^{π/2} F(q,α)^2 sinα dα.

    Данные: amp_fn, q (нм^-1), V
    Проверка: нет независимой проверки (вспомогательная)
    Используется в: src/FormFactors/form_factors.py::cylinder_intensity, ::core_shell_cylinder_intensity, ::cylinder_radial_intensity
    """
    out = np.empty_like(np.asarray(q, dtype=float))
    for i, qi in enumerate(np.asarray(q, dtype=float)):
        val, _ = quad(lambda a: amp_fn(qi, a) ** 2 * np.sin(a), 0.0, np.pi / 2.0, limit=500)
        out[i] = val / V
    return out


# ---------------------------------------------------------------------------
# Сферы
# ---------------------------------------------------------------------------

def sphere_amplitude(q: np.ndarray, R: float, rho: float = 1.0) -> np.ndarray:
    """Амплитуда однородной сферы радиуса R с контрастом rho.

    Данные: q (нм^-1), R, rho
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: src/FormFactors/form_factors.py::sphere_intensity
    """
    V = 4.0 / 3.0 * np.pi * R ** 3
    return rho * V * _sphere_f(np.asarray(q, dtype=float) * R)


def sphere_intensity(q: np.ndarray, R: float, rho: float = 1.0) -> np.ndarray:
    """P(q) однородной сферы, нормированная на объём.

    Данные: q (нм^-1), R, rho
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: src/FormFactors/form_factors.py::sphere_intensity_pd
    """
    V = 4.0 / 3.0 * np.pi * R ** 3
    return sphere_amplitude(q, R, rho) ** 2 / V


def core_shell_sphere_amplitude(
    q: np.ndarray,
    R_core: float,
    R_shell: float,
    rho_core: float,
    rho_shell: float,
    rho_solv: float = 0.0,
) -> np.ndarray:
    """Амплитуда сферы с ядром и оболочкой (концентрические сферы).

    Данные: q (нм^-1), R_core, R_shell, rho_core, rho_shell, rho_solv
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: src/FormFactors/form_factors.py::core_shell_sphere_intensity
    """
    Vc = 4.0 / 3.0 * np.pi * R_core ** 3
    Vs = 4.0 / 3.0 * np.pi * R_shell ** 3
    q = np.asarray(q, dtype=float)
    return (
        (rho_core - rho_shell) * Vc * _sphere_f(q * R_core)
        + (rho_shell - rho_solv) * Vs * _sphere_f(q * R_shell)
    )


def core_shell_sphere_intensity(
    q: np.ndarray,
    R_core: float,
    R_shell: float,
    rho_core: float,
    rho_shell: float,
    rho_solv: float = 0.0,
) -> np.ndarray:
    """P(q) сферы с ядром и оболочкой, нормированная на объём оболочки.

    Данные: q (нм^-1), R_core, R_shell, rho_core, rho_shell, rho_solv
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: экспортируется в src/FormFactors/__init__.py
    """
    Vs = 4.0 / 3.0 * np.pi * R_shell ** 3
    return core_shell_sphere_amplitude(q, R_core, R_shell, rho_core, rho_shell, rho_solv) ** 2 / Vs


# ---------------------------------------------------------------------------
# Цилиндры
# ---------------------------------------------------------------------------

def infinite_cylinder_amplitude(q: np.ndarray, R: float, rho: float = 1.0, L: float = 1.0) -> np.ndarray:
    """Амплитуда бесконечного цилиндра (L >> R), ось ⊥ q.

    F(q) = 2 V rho J1(qR)/(qR), V = π R^2 L.

    Данные: q (нм^-1), R, rho, L
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: src/FormFactors/form_factors.py::infinite_cylinder_intensity
    """
    V = np.pi * R ** 2 * L
    return 2.0 * rho * V * _j1_over_x(np.asarray(q, dtype=float) * R)


def infinite_cylinder_intensity(q: np.ndarray, R: float, rho: float = 1.0, L: float = 1.0) -> np.ndarray:
    """P(q) бесконечного цилиндра, нормированная на объём.

    Ориентационно-усреднённый форм-фактор длинного цилиндра (L >> R):
    P(q) = rho^2 V [2 J1(qR)/(qR)]^2 / q,  V = pi R^2 L.
    Фактор 1/q — от усреднения по ориентации оси цилиндра.

    Данные: q (нм^-1), R, rho, L
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: src/FormFactors/form_factors.py::infinite_cylinder_intensity_pd
    """
    V = np.pi * R ** 2 * L
    q = np.asarray(q, dtype=float)
    return infinite_cylinder_amplitude(q, R, rho, L) ** 2 / V / q


def cylinder_amplitude(
    q: np.ndarray,
    R: float,
    L: float,
    rho: float = 1.0,
    alpha: float | None = None,
) -> np.ndarray | Callable[[float, float], float]:
    """Амплитуда конечного цилиндра.

    alpha — угол между осью цилиндра и q. Если alpha=None, возвращается
    функция F(q, alpha) для последующего усреднения по ориентации.

    Данные: q (нм^-1), R, L, rho, alpha
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: src/FormFactors/form_factors.py::cylinder_intensity
    """
    V = np.pi * R ** 2 * L

    def amp(qi: float, a: float) -> float:
        x = qi * np.sin(a) * R
        y = qi * np.cos(a) * L / 2.0
        return 2.0 * rho * V * _j1_over_x(x) * _sinc(y)

    if alpha is None:
        return amp
    return amp(np.asarray(q, dtype=float), alpha)


def cylinder_intensity(q: np.ndarray, R: float, L: float, rho: float = 1.0) -> np.ndarray:
    """P(q) конечного цилиндра, усреднённая по ориентации.

    Данные: q (нм^-1), R, L, rho
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: экспортируется в src/FormFactors/__init__.py
    """
    V = np.pi * R ** 2 * L
    return _orient_avg(cylinder_amplitude(q, R, L, rho), q, V)


def core_shell_cylinder_amplitude(
    q: np.ndarray,
    R_core: float,
    R_shell: float,
    L: float,
    rho_core: float,
    rho_shell: float,
    rho_solv: float = 0.0,
    alpha: float | None = None,
) -> np.ndarray | Callable[[float, float], float]:
    """Амплитуда цилиндра с ядром и оболочкой (одинаковая длина L).

    Данные: q (нм^-1), R_core, R_shell, L, rho_core, rho_shell, rho_solv, alpha
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: src/FormFactors/form_factors.py::core_shell_cylinder_intensity
    """
    Vc = np.pi * R_core ** 2 * L
    Vs = np.pi * R_shell ** 2 * L

    def amp(qi: float, a: float) -> float:
        xc = qi * np.sin(a) * R_core
        xs = qi * np.sin(a) * R_shell
        y = qi * np.cos(a) * L / 2.0
        return (
            2.0 * (rho_core - rho_shell) * Vc * _j1_over_x(xc) * _sinc(y)
            + 2.0 * (rho_shell - rho_solv) * Vs * _j1_over_x(xs) * _sinc(y)
        )

    if alpha is None:
        return amp
    return amp(np.asarray(q, dtype=float), alpha)


def core_shell_cylinder_intensity(
    q: np.ndarray,
    R_core: float,
    R_shell: float,
    L: float,
    rho_core: float,
    rho_shell: float,
    rho_solv: float = 0.0,
) -> np.ndarray:
    """P(q) цилиндра с ядром и оболочкой, усреднённая по ориентации.

    Данные: q (нм^-1), R_core, R_shell, L, rho_core, rho_shell, rho_solv
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: экспортируется в src/FormFactors/__init__.py
    """
    Vs = np.pi * R_shell ** 2 * L
    return _orient_avg(
        core_shell_cylinder_amplitude(q, R_core, R_shell, L, rho_core, rho_shell, rho_solv),
        q, Vs,
    )


# ---------------------------------------------------------------------------
# Цилиндр с произвольным радиальным профилем электронной плотности
# ---------------------------------------------------------------------------

def cylinder_radial_amplitude(
    q: np.ndarray,
    rho_r: Callable[[float], float],
    R: float,
    L: float,
    alpha: float | None = None,
) -> np.ndarray | Callable[[float, float], float]:
    """Форм-фактор цилиндра с радиальным профилем плотности rho_r(r).

    F(q,α) = 2π L sinc(qL cosα/2) ∫_0^R rho_r(r) J0(q sinα r) r dr

    rho_r — callable: r -> электронная плотность (e-/nm^3).
    alpha=None → возвращает функцию F(q, alpha) для усреднения по ориентации.

    Данные: q (нм^-1), rho_r, R, L, alpha
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: src/FormFactors/form_factors.py::cylinder_radial_intensity
    """
    def amp(qi: float, a: float) -> float:
        s = np.sin(a)
        integ, _ = quad(lambda r: rho_r(r) * j0(qi * s * r) * r, 0.0, R, limit=200)
        return 2.0 * np.pi * L * _sinc(qi * np.cos(a) * L / 2.0) * integ

    if alpha is None:
        return amp
    return amp(np.asarray(q, dtype=float), alpha)


def cylinder_radial_intensity(
    q: np.ndarray,
    rho_r: Callable[[float], float],
    R: float,
    L: float,
) -> np.ndarray:
    """P(q) цилиндра с радиальным профилем, усреднённая по ориентации.

    Нормировка на объём цилиндра V = π R^2 L.

    Данные: q (нм^-1), rho_r, R, L
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: экспортируется в src/FormFactors/__init__.py
    """
    V = np.pi * R ** 2 * L
    return _orient_avg(cylinder_radial_amplitude(q, rho_r, R, L), q, V)


def infinite_cylinder_radial_amplitude(
    q: np.ndarray,
    rho_r: Callable[[float], float],
    R: float,
    L: float = 1.0,
) -> np.ndarray:
    """Форм-фактор бесконечного цилиндра с радиальным профилем rho_r(r).

    F(q) = 2π L ∫_0^R rho_r(r) J0(q r) r dr  (ось ⊥ q).

    Данные: q (нм^-1), rho_r, R, L
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: src/FormFactors/form_factors.py::infinite_cylinder_radial_intensity
    """
    q = np.asarray(q, dtype=float)
    out = np.empty_like(q)
    for i, qi in enumerate(q):
        integ, _ = quad(lambda r: rho_r(r) * j0(qi * r) * r, 0.0, R, limit=200)
        out[i] = 2.0 * np.pi * L * integ
    return out


def infinite_cylinder_radial_intensity(
    q: np.ndarray,
    rho_r: Callable[[float], float],
    R: float,
    L: float = 1.0,
) -> np.ndarray:
    """P(q) бесконечного цилиндра с радиальным профилем, нормированная на объём.

    Данные: q (нм^-1), rho_r, R, L
    Проверка: задача 2 (фит форм-факторов) — checks/check.py, погрешность ≤ 5%
    Используется в: экспортируется в src/FormFactors/__init__.py
    """
    V = np.pi * R ** 2 * L
    return infinite_cylinder_radial_amplitude(q, rho_r, R, L) ** 2 / V

# ---------------------------------------------------------------------------
# Полидисперсность радиуса (Gaussian)
# ---------------------------------------------------------------------------

def _gauss_pd_weight(r: np.ndarray, R0: float, sigma: float) -> np.ndarray:
    """Вес Gaussian-распределения радиуса: exp(-(r-R0)^2/(2 sigma^2)).

    sigma = pd * R0, где pd — относительная полидисперсность.
    """
    return np.exp(-0.5 * ((r - R0) / sigma) ** 2)


def sphere_intensity_pd(
    q: np.ndarray,
    R0: float,
    pd: float = 0.0,
    rho: float = 1.0,
    n_pts: int = 200,
    r_lo: float = 0.0,
    r_hi: float | None = None,
) -> np.ndarray:
    """P(q) сферы с Gaussian-полидисперсностью радиуса.

    Усреднение по радиусу: <P(q)> = ∫ P(q,R) w(R) dR / ∫ w(R) dR,
    где w(R) — гауссово распределение со средним R0 и sigma = pd*R0.
    """
    sigma = pd * R0
    if sigma <= 0:
        return sphere_intensity(q, R0, rho)
    if r_hi is None:
        r_hi = R0 + 4.0 * sigma
    r = np.linspace(r_lo, r_hi, n_pts)
    w = _gauss_pd_weight(r, R0, sigma)
    # P(q,R) = rho^2 V(R) f(qR)^2  (V = 4/3 pi R^3)
    V = 4.0 / 3.0 * np.pi * r ** 3
    f = _sphere_f(np.asarray(q, dtype=float)[:, None] * r[None, :])
    P = rho ** 2 * V[None, :] * f ** 2
    return np.trapezoid(P * w[None, :], r, axis=1) / np.trapezoid(w, r)


def infinite_cylinder_intensity_pd(
    q: np.ndarray,
    R0: float,
    pd: float = 0.0,
    rho: float = 1.0,
    L: float = 1.0,
    n_pts: int = 200,
    r_lo: float = 0.0,
    r_hi: float | None = None,
) -> np.ndarray:
    """P(q) бесконечного цилиндра с Gaussian-полидисперсностью радиуса.

    Усреднение по радиусу: <P(q)> = ∫ P(q,R) w(R) dR / ∫ w(R) dR.
    P(q,R) = rho^2 V(R) [2 J1(qR)/(qR)]^2 / q, V = pi R^2 L.

    Фактор 1/q — от ориентационного усреднения длинного цилиндра (L >> R).
    """
    sigma = pd * R0
    if sigma <= 0:
        return infinite_cylinder_intensity(q, R0, rho, L)
    if r_hi is None:
        r_hi = R0 + 4.0 * sigma
    r = np.linspace(r_lo, r_hi, n_pts)
    w = _gauss_pd_weight(r, R0, sigma)
    V = np.pi * r ** 2 * L
    j = 2.0 * _j1_over_x(np.asarray(q, dtype=float)[:, None] * r[None, :])
    P = rho ** 2 * V[None, :] * j ** 2 / np.asarray(q, dtype=float)[:, None]
    return np.trapezoid(P * w[None, :], r, axis=1) / np.trapezoid(w, r)


# ---------------------------------------------------------------------------
# Нормированные форм-факторы (P(0) = 1) с полидисперсностью
# ---------------------------------------------------------------------------

def sphere_intensity_pd_norm(
    q: np.ndarray,
    R0: float,
    pd: float = 0.0,
    n_pts: int = 200,
    r_lo: float = 0.0,
    r_hi: float | None = None,
) -> np.ndarray:
    """P(q) сферы с Gaussian-полидисперсностью, нормированный на P(0)=1.

    <P(q)> = ∫ f(qR)^2 w(R) dR / ∫ w(R) dR,  f(x)=3(sin x - x cos x)/x^3.
    """
    sigma = pd * R0
    if sigma <= 0:
        return _sphere_f(np.asarray(q, dtype=float) * R0) ** 2
    if r_hi is None:
        r_hi = R0 + 4.0 * sigma
    r = np.linspace(r_lo, r_hi, n_pts)
    w = _gauss_pd_weight(r, R0, sigma)
    f = _sphere_f(np.asarray(q, dtype=float)[:, None] * r[None, :])
    return np.trapezoid(f ** 2 * w[None, :], r, axis=1) / np.trapezoid(w, r)


def infinite_cylinder_intensity_pd_norm(
    q: np.ndarray,
    R0: float,
    pd: float = 0.0,
    n_pts: int = 200,
    r_lo: float = 0.0,
    r_hi: float | None = None,
) -> np.ndarray:
    """P(q) бесконечного цилиндра с полидисперсностью, нормированный на P(0)=1.

    Ориентационно-усреднённый форм-фактор длинного цилиндра (L >> R) содержит
    фактор 1/q от усреднения по углу между осью цилиндра и q:

        <P(q)> = (1/q) ∫ [2 J1(qR)/(qR)]^2 w(R) dR / ∫ w(R) dR.

    Без фактора 1/q асимптотика была бы [2J1(qR)/(qR)]^2 ~ q^-3, тогда как
    правильная асимптотика длинного цилиндра ~ q^-4. Фактор 1/q критичен для
    корректного восстановления радиуса при фитировании.
    """
    q = np.asarray(q, dtype=float)
    sigma = pd * R0
    if sigma <= 0:
        return (2.0 * _j1_over_x(q * R0)) ** 2 / q
    if r_hi is None:
        r_hi = R0 + 4.0 * sigma
    r = np.linspace(r_lo, r_hi, n_pts)
    w = _gauss_pd_weight(r, R0, sigma)
    j = 2.0 * _j1_over_x(q[:, None] * r[None, :])
    return np.trapezoid(j ** 2 * w[None, :], r, axis=1) / np.trapezoid(w, r) / q
