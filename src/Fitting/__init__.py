"""Фитирование SAXS-кривых форм-факторами P(q)."""
from .fit_form_factors import (
    FitResult,
    fit_form_factor,
    fit_from_file,
    fit_result_to_dict,
)

__all__ = [
    "FitResult",
    "fit_form_factor",
    "fit_from_file",
    "fit_result_to_dict",
]