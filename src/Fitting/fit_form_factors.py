"""Фитирование SAXS-кривых форм-факторами P(q).

Модель: I(q) = scale * P(q; R, pd) + background.

Форм-факторы с Gaussian-полидисперсностью берутся из
`src.FormFactors.form_factors`.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.optimize import least_squares

from FormFactors.form_factors import (
    sphere_intensity_pd_norm,
    infinite_cylinder_intensity_pd_norm,
)


@dataclass
class FitResult:
    """Результат фита одного форм-фактора."""

    model: str
    radius_nm: float
    polydispersity: float
    scale: float
    background: float
    q_nm: np.ndarray
    I_exp: np.ndarray
    I_fit: np.ndarray
    success: bool
    message: str


def _model_intensity(
    q_nm: np.ndarray,
    radius_nm: float,
    scale: float,
    background: float,
    pd: float,
    model: str,
) -> np.ndarray:
    """Рассчитать I(q) = scale * P(q) + background для заданной модели."""
    if model == "sphere":
        P = sphere_intensity_pd_norm(q_nm, radius_nm, pd)
    elif model == "cylinder":
        P = infinite_cylinder_intensity_pd_norm(q_nm, radius_nm, pd)
    else:
        raise ValueError(f"unknown model: {model}")
    return scale * P + background


def _residual(
    params: np.ndarray,
    q_nm: np.ndarray,
    I_exp: np.ndarray,
    pd: float,
    model: str,
) -> np.ndarray:
    """Невязка в лог-пространстве (относительная), устойчива к масштабу."""
    radius_nm, scale, background = params
    I_fit = _model_intensity(q_nm, radius_nm, scale, background, pd, model)
    # относительная невязка в лог-пространстве
    return np.log(I_fit + 1e-12) - np.log(I_exp + 1e-12)


def _residual_fixed_bg(
    params: np.ndarray,
    q_nm: np.ndarray,
    I_exp: np.ndarray,
    pd: float,
    model: str,
    background: float,
) -> np.ndarray:
    """Невязка в лог-пространстве с фиксированным фоном: (radius_nm, scale)."""
    radius_nm, scale = params
    I_fit = _model_intensity(q_nm, radius_nm, scale, background, pd, model)
    return np.log(I_fit + 1e-12) - np.log(I_exp + 1e-12)


def fit_form_factor(
    q_ang: np.ndarray,
    I_exp: np.ndarray,
    model: str,
    pd: float = 0.1,
    radius_guess: float | None = None,
    scale_guess: float | None = None,
    background_guess: float = 0.0,
    bounds: tuple[tuple[float, float, float], tuple[float, float, float]] | None = None,
    q_range: tuple[float, float] | None = None,
    fix_background: float | None = None,
) -> FitResult:
    """Фитировать форм-фактор по экспериментальной кривой.

    Параметры:
        q_ang: q в ангстремах^-1 (как в SasView).
        I_exp: интенсивность.
        model: "sphere" или "cylinder".
        pd: фиксированная полидисперсность радиуса.
        radius_guess: начальное приближение радиуса (нм). Если None —
            оценивается по первому минимуму кривой.
        bounds: (lower, upper) для (radius_nm, scale, background).
        q_range: (q_min, q_max) в нм^-1 — ограничить фит областью осцилляций.
            Если None — фитируется весь диапазон.
        fix_background: если задано, фон фиксируется на этом значении
            (не фитируется). Полезно, когда фон оценивается по хвосту кривой.

    Возвращает FitResult с параметрами и рассчитанной кривой.
    """
    q_nm = np.asarray(q_ang, dtype=float) * 10.0
    I_exp = np.asarray(I_exp, dtype=float)

    # ограничить область фита
    if q_range is not None:
        qmin, qmax = q_range
        mask = (q_nm >= qmin) & (q_nm <= qmax)
        q_fit = q_nm[mask]
        I_fit_exp = I_exp[mask]
    else:
        q_fit = q_nm
        I_fit_exp = I_exp

    # Оценка начального радиуса по первому минимуму.
    # Сфера: первый минимум f(qR) при qR ~ 4.493.
    # Длинный цилиндр: первый минимум [2J1(qR)/(qR)]^2 при qR ~ 3.832
    #   (первый ноль J1).
    if radius_guess is None:
        # ищем первый локальный минимум лог-интенсивности
        logI = np.log(I_fit_exp + 1e-12)
        d = np.diff(logI)
        sign = np.sign(d)
        # первый переход с убывания на возрастание
        idx = np.where((sign[:-1] < 0) & (sign[1:] > 0))[0]
        if len(idx) > 0:
            q_min = q_fit[idx[0] + 1]
            if model == "cylinder":
                radius_guess = 3.8317 / q_min if q_min > 0 else 1.0
            else:
                radius_guess = 4.493 / q_min if q_min > 0 else 1.0
        else:
            radius_guess = 1.0

    if scale_guess is None:
        scale_guess = I_fit_exp[0] / (sphere_intensity_pd_norm(q_fit[:1], radius_guess, pd)[0] + 1e-12)

    if fix_background is not None:
        # фон фиксирован — фитируем только радиус и scale
        bg = float(fix_background)
        if bounds is None:
            bounds = (
                (radius_guess * 0.1, scale_guess * 1e-3),
                (radius_guess * 10.0, scale_guess * 1e3),
            )
        p0 = np.array([radius_guess, scale_guess], dtype=float)
        lo = np.array(bounds[0], dtype=float)
        hi = np.array(bounds[1], dtype=float)
        res = least_squares(
            _residual_fixed_bg,
            p0,
            args=(q_fit, I_fit_exp, pd, model, bg),
            bounds=(lo, hi),
            max_nfev=2000,
        )
        radius_nm, scale = res.x
        I_fit = _model_intensity(q_nm, radius_nm, scale, bg, pd, model)
        return FitResult(
            model=model,
            radius_nm=float(radius_nm),
            polydispersity=float(pd),
            scale=float(scale),
            background=float(bg),
            q_nm=q_nm,
            I_exp=I_exp,
            I_fit=I_fit,
            success=res.success,
            message=res.message,
        )

    if bounds is None:
        bounds = (
            (radius_guess * 0.1, scale_guess * 1e-3, -1e-3),
            (radius_guess * 10.0, scale_guess * 1e3, I_fit_exp.max()),
        )

    p0 = np.array([radius_guess, scale_guess, background_guess], dtype=float)
    lo = np.array(bounds[0], dtype=float)
    hi = np.array(bounds[1], dtype=float)

    res = least_squares(
        _residual,
        p0,
        args=(q_fit, I_fit_exp, pd, model),
        bounds=(lo, hi),
        max_nfev=2000,
    )

    radius_nm, scale, background = res.x
    I_fit = _model_intensity(q_nm, radius_nm, scale, background, pd, model)

    return FitResult(
        model=model,
        radius_nm=float(radius_nm),
        polydispersity=float(pd),
        scale=float(scale),
        background=float(background),
        q_nm=q_nm,
        I_exp=I_exp,
        I_fit=I_fit,
        success=res.success,
        message=res.message,
    )


def fit_from_file(
    data_path: Path,
    model: str,
    pd: float = 0.1,
    radius_guess: float | None = None,
    q_range: tuple[float, float] | None = None,
) -> FitResult:
    """Фитировать кривую из текстового файла (колонки q, I)."""
    data = np.loadtxt(data_path, comments="#", skiprows=1)
    q_ang = data[:, 0]
    I_exp = data[:, 1]
    return fit_form_factor(q_ang, I_exp, model, pd=pd, radius_guess=radius_guess, q_range=q_range)


def fit_result_to_dict(result: FitResult) -> dict:
    """Сериализовать результат в dict для JSON (без массивов)."""
    return {
        "model": result.model,
        "parameters": {
            "radius_nm": result.radius_nm,
            "polydispersity": result.polydispersity,
            "scale": result.scale,
            "background": result.background,
        },
        "success": result.success,
        "message": result.message,
    }