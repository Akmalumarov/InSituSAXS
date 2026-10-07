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
    """sin(x)/x, устойчиво при x -> 0."""
    return np.sinc(np.asarray(x, dtype=float) / np.pi)


def _j1_over_x(x: np.ndarray) -> np.ndarray:
    """J1(x)/x, устойчиво при x -> 0 (предел = 1/2)."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    out[small] = 0.5 - x[small] ** 2 / 16.0  # разложение J1(x)/x
    out[~small] = j1(x[~small]) / x[~small]
    return out


def _sphere_f(x: np.ndarray) -> np.ndarray:
    """f(x) = 3 (sin x - x cos x) / x^3, f(0) = 1."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    out[small] = 1.0 - x[small] ** 2 / 10.0  # разложение
    out[~small] = 3.0 * (np.sin(x[~small]) - x[~small] * np.cos(x[~small])) / x[~small] ** 3
    return out


def _orient_avg(amp_fn: Callable[[float], float], q: np.ndarray, V: float) -> np.ndarray:
    """Усреднение по ориентации: P(q) = (1/V) ∫_0^{π/2} F(q,α)^2 sinα dα."""
    out = np.empty_like(np.asarray(q, dtype=float))
    for i, qi in enumerate(np.asarray(q, dtype=float)):
        val, _ = quad(lambda a: amp_fn(qi, a) ** 2 * np.sin(a), 0.0, np.pi / 2.0, limit=500)
        out[i] = val / V
    return out


# ---------------------------------------------------------------------------
# Сферы
# ---------------------------------------------------------------------------

def sphere_amplitude(q: np.ndarray, R: float, rho: float = 1.0) -> np.ndarray:
    """Амплитуда однородной сферы радиуса R с контрастом rho."""
    V = 4.0 / 3.0 * np.pi * R ** 3
    return rho * V * _sphere_f(np.asarray(q, dtype=float) * R)


def sphere_intensity(q: np.ndarray, R: float, rho: float = 1.0) -> np.ndarray:
    """P(q) однородной сферы, нормированная на объём."""
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
    """Амплитуда сферы с ядром и оболочкой (концентрические сферы)."""
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
    """P(q) сферы с ядром и оболочкой, нормированная на объём оболочки."""
    Vs = 4.0 / 3.0 * np.pi * R_shell ** 3
    return core_shell_sphere_amplitude(q, R_core, R_shell, rho_core, rho_shell, rho_solv) ** 2 / Vs


# ---------------------------------------------------------------------------
# Цилиндры
# ---------------------------------------------------------------------------

def infinite_cylinder_amplitude(q: np.ndarray, R: float, rho: float = 1.0, L: float = 1.0) -> np.ndarray:
    """Амплитуда бесконечного цилиндра (L >> R), ось ⊥ q.

    F(q) = 2 V rho J1(qR)/(qR), V = π R^2 L.
    """
    V = np.pi * R ** 2 * L
    return 2.0 * rho * V * _j1_over_x(np.asarray(q, dtype=float) * R)


def infinite_cylinder_intensity(q: np.ndarray, R: float, rho: float = 1.0, L: float = 1.0) -> np.ndarray:
    """P(q) бесконечного цилиндра, нормированная на объём."""
    V = np.pi * R ** 2 * L
    return infinite_cylinder_amplitude(q, R, rho, L) ** 2 / V


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
    """P(q) конечного цилиндра, усреднённая по ориентации."""
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
    """Амплитуда цилиндра с ядром и оболочкой (одинаковая длина L)."""
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
    """P(q) цилиндра с ядром и оболочкой, усреднённая по ориентации."""
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
    """P(q) бесконечного цилиндра с радиальным профилем, нормированная на объём."""
    V = np.pi * R ** 2 * L
    return infinite_cylinder_radial_amplitude(q, rho_r, R, L) ** 2 / V