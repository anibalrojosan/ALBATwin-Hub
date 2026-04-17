"""
Algebraic kinetic modifiers f_T, f_pH, f_I, f_DO (``docs/MATH_MODEL.md`` §3).

Pure functions; no model state beyond scalar arguments.
"""

from __future__ import annotations

from bioprocess_twin.models.kinetic_parameters import CardinalPH, CardinalTemperature

_DENOM_EPS = 1e-18


def f_t_growth_ctmi(t_c: float, card: CardinalTemperature) -> float:
    """
    CTMI growth factor (Rosso et al., 1993), ``MATH_MODEL.md`` §3.1.

    Returns 0 outside ``[t_min, t_max]`` (inclusive bounds yield 0 in numerator).
    """
    t_min = card.t_min
    t_opt = card.t_opt
    t_max = card.t_max
    if t_c < t_min or t_c > t_max:
        return 0.0
    num = (t_c - t_max) * (t_c - t_min) ** 2
    inner = (t_opt - t_min) * (t_c - t_opt) - (t_opt - t_max) * (t_opt + t_min - 2.0 * t_c)
    den = (t_opt - t_min) * inner
    if abs(den) < _DENOM_EPS:
        return 0.0
    return float(num / den)


def f_t_decay_arrhenius(t_c: float, theta: float) -> float:
    """``θ^(T−20)`` decay temperature factor, §3.1."""
    if theta <= 0.0:
        raise ValueError("theta must be positive for real-valued decay factor")
    return float(theta ** (t_c - 20.0))


def f_ph_cpm(ph: float, card: CardinalPH) -> float:
    """
    Cardinal pH model (CPM), §3.2.

    Returns 0 if ``ph`` is outside ``[ph_min, ph_max]``.
    """
    lo, hi = card.ph_min, card.ph_max
    if ph < lo or ph > hi:
        return 0.0
    num = (ph - lo) * (ph - hi)
    den = num - (ph - card.ph_opt) ** 2
    if abs(den) < _DENOM_EPS:
        return 0.0
    return float(num / den)


def f_i_haldane(i: float, mu_max: float, alpha: float, i_opt: float) -> float:
    """
    Haldane-type light response (Bernard & Rémond, 2012), §3.3.

    At ``I == i_opt`` the factor equals ``mu_max``.
    """
    if i < 0.0:
        raise ValueError("irradiance I must be non-negative")
    if mu_max <= 0.0 or alpha <= 0.0 or i_opt <= 0.0:
        raise ValueError("mu_max, alpha, and i_opt must be positive")
    r = i / i_opt - 1.0
    return float(mu_max / (1.0 + (mu_max / alpha) * r * r))


def f_do_growth_hill(s_o2: float, ec50: float, n: float) -> float:
    """Algal growth oxygen inhibition (Hill), §3.4 first equation."""
    if s_o2 < 0.0 or ec50 <= 0.0 or n <= 0.0:
        raise ValueError("s_o2 must be >= 0; ec50 and n must be positive")
    return float(ec50**n / (s_o2**n + ec50**n))


def f_do_decay_hill(s_o2: float, ec50: float, n: float) -> float:
    """Algal decay oxygen enhancement (Hill), §3.4 second equation."""
    if s_o2 < 0.0 or ec50 <= 0.0 or n <= 0.0:
        raise ValueError("s_o2 must be >= 0; ec50 and n must be positive")
    return float(s_o2**n / (s_o2**n + ec50**n))
