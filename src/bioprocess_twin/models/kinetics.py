"""
ALBA biological process rates ρ.

Process order matches ``get_petersen_matrix()`` row index ``i`` = ρ_{i+1}.
Expressions: ``docs/MATH_MODEL.md`` §4; modifiers §3.

§4 writes ρ₁ = μ_max,ALG · f_I · … with f_I from §3.3; §3.3 defines f_I such that
f_I(I_opt) = μ_max. To avoid μ_max² at optimum, growth uses
``f_i_haldane(...) / μ_max,ALG`` as a dimensionless light factor (see docstring below).
"""

from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from bioprocess_twin.core.state import StateVector, StateVectorVariant

from .kinetic_modifiers import (
    f_do_decay_hill,
    f_do_growth_hill,
    f_i_haldane,
    f_ph_cpm,
    f_t_decay_arrhenius,
    f_t_growth_ctmi,
)
from .kinetic_parameters import KineticParameters, default_alba
from .stoichiometry import (
    N_PROCESSES,
    N_STATE,
    S_IC,
    S_ND,
    S_NH,
    S_NO2,
    S_NO3,
    S_O2,
    S_PO4,
    S_S,
    X_ALG,
    X_AOB,
    X_H,
    X_NOB,
    X_S,
)


class EnvConditions(BaseModel):
    """Environmental inputs for kinetic modifiers (§3) and process rates (§4)."""

    model_config = ConfigDict(frozen=True)

    temperature_C: float = Field(..., description="Temperature [°C]; may be negative (CTMI T_min for algae).")
    pH: float = Field(..., description="pH for cardinal pH factor (CPM).")
    irradiance_umol_m2_s: float = Field(
        ...,
        ge=0.0,
        description="Photosynthetically active irradiance [μmol m⁻² s⁻¹] for algal f_I (§3.3).",
    )


def _nonneg(x: float) -> float:
    return float(x) if x > 0.0 else 0.0


def _monod(s: float, k: float) -> float:
    """s / (k + s) with s clamped non-negative; k > 0 assumed."""
    s = _nonneg(s)
    return s / (k + s)


def _anox_switch(s_o2: float, k_o: float) -> float:
    """K_O / (K_O + S_O2) for anoxic processes."""
    return k_o / (k_o + _nonneg(s_o2))


def calculate_rates(
    state: np.ndarray | StateVector,
    env_conditions: EnvConditions,
    kinetic_parameters: KineticParameters | None = None,
) -> np.ndarray:
    """
    Return the vector of **19** process rates ρ₁…ρ₁₉ as in ``MATH_MODEL.md`` §4.

    Each ρᵢ has the form [d⁻¹] × (concentration products) so that ``Sᵀρ`` matches
    the SI mass-rate convention.

    **Light factor:** §3.3 defines f_I(I) with f_I(I_opt) = μ_max; §4 multiplies
    μ_max,ALG · f_I. The implementation uses
    ``f_i_dim = f_i_haldane(I, μ_max,ALG, α, I_opt) / μ_max,ALG`` in ρ₁ and ρ₂ so
    that f_i_dim(I_opt) = 1.

    State layout: first **17** entries are Casagli SI variables. A length-**18**
    ndarray uses only ``state[:17]`` (proton closure is ignored for ρ).

    Args:
        state: ``(17,)`` or ``(18,)`` concentrations, or ``StateVector``.
        env_conditions: Temperature, pH, irradiance.
        kinetic_parameters: Defaults to ``default_alba()`` when omitted.

    Returns:
        ``ndarray`` shape ``(19,)``, ``dtype`` float64.

    Raises:
        ValueError: If ``state`` array size is not 17 or 18.
    """
    params = kinetic_parameters if kinetic_parameters is not None else default_alba()

    if isinstance(state, StateVector):
        conc = state.to_array(variant=StateVectorVariant.SI)
    else:
        flat = np.asarray(state, dtype=np.float64).ravel()
        n = flat.size
        if n == N_STATE:
            conc = flat
        elif n == N_STATE + 1:
            conc = flat[:N_STATE]
        else:
            raise ValueError(f"state as ndarray must have length {N_STATE} or {N_STATE + 1} (proton layout), got {n}")

    T = env_conditions.temperature_C
    ph = env_conditions.pH
    irr = env_conditions.irradiance_umol_m2_s

    # --- Modifiers by group (§3, §4 table) ---
    f_t_alg = f_t_growth_ctmi(T, params.cardinal_temp_alg)
    f_ph_alg = f_ph_cpm(ph, params.cardinal_ph_alg)
    f_t_h = f_t_growth_ctmi(T, params.cardinal_temp_h)
    f_ph_h = f_ph_cpm(ph, params.cardinal_ph_h)
    f_t_aob = f_t_growth_ctmi(T, params.cardinal_temp_aob)
    f_ph_aob = f_ph_cpm(ph, params.cardinal_ph_aob)
    f_t_nob = f_t_growth_ctmi(T, params.cardinal_temp_nob)
    f_ph_nob = f_ph_cpm(ph, params.cardinal_ph_nob)

    f_i_dim = f_i_haldane(irr, params.mu_max_alg, params.alpha_light, params.i_opt) / params.mu_max_alg
    s_o2 = float(conc[S_O2])
    f_do_g = f_do_growth_hill(s_o2, params.ec50_o2, params.hill_n_o2)
    f_do_d = f_do_decay_hill(s_o2, params.ec50_o2, params.hill_n_o2)

    theta_alg_T = f_t_decay_arrhenius(T, params.theta_alg)
    theta_h_T = f_t_decay_arrhenius(T, params.theta_h)
    theta_hyd_T = f_t_decay_arrhenius(T, params.theta_hyd)
    theta_amm_T = f_t_decay_arrhenius(T, params.theta_amm)
    theta_aob_T = f_t_decay_arrhenius(T, params.theta_aob)
    theta_nob_T = f_t_decay_arrhenius(T, params.theta_nob)

    x_alg = float(conc[X_ALG])
    x_aob = float(conc[X_AOB])
    x_nob = float(conc[X_NOB])
    x_h = float(conc[X_H])
    x_s = float(conc[X_S])

    rho = np.zeros(N_PROCESSES, dtype=np.float64)

    # --- Algae ρ1–ρ4 (§4); K_C,K_N,K_NO3,K_P,K_O → k_c_alg, k_n_alg, … ---
    m_ic_alg = _monod(float(conc[S_IC]), params.k_c_alg)
    m_nh_alg = _monod(float(conc[S_NH]), params.k_n_alg)
    m_no3_alg = _monod(float(conc[S_NO3]), params.k_no3_alg)
    m_po4_alg = _monod(float(conc[S_PO4]), params.k_p_alg)
    m_o2_alg = _monod(s_o2, params.k_o_alg)

    liebig_rho1 = min(m_ic_alg, m_nh_alg, m_po4_alg)
    rho[0] = params.mu_max_alg * f_i_dim * f_t_alg * f_ph_alg * f_do_g * liebig_rho1 * x_alg

    kn_knh_alg = _monod(float(conc[S_NH]), params.k_n_alg)  # S_NH/(K_N+S_NH)
    rho[1] = (
        params.mu_max_alg
        * f_i_dim
        * f_t_alg
        * f_ph_alg
        * f_do_g
        * (1.0 - kn_knh_alg)  # K_N/(K_N+S_NH) = 1 - S_NH/(K_N+S_NH)
        * min(m_ic_alg, m_no3_alg, m_po4_alg)
        * x_alg
    )

    rho[2] = params.b_max_r_alg * f_t_alg * f_ph_alg * m_o2_alg * x_alg
    rho[3] = params.b_max_d_alg * theta_alg_T * f_ph_alg * f_do_d * x_alg

    # --- Heterotrophs ρ5–ρ10 ---
    m_ss = _monod(float(conc[S_S]), params.k_s_h)
    m_o2_h = _monod(s_o2, params.k_o_h)
    m_nh_h = _monod(float(conc[S_NH]), params.k_n_h)
    m_no2_h = _monod(float(conc[S_NO2]), params.k_no2_h)
    m_no3_h = _monod(float(conc[S_NO3]), params.k_no3_h)
    m_po4_h = _monod(float(conc[S_PO4]), params.k_p_h)
    kn_knh_h = _monod(float(conc[S_NH]), params.k_n_h)
    anox = _anox_switch(s_o2, params.k_o_h)

    rho[4] = params.mu_max_h * f_t_h * f_ph_h * min(m_ss, m_o2_h, m_nh_h, m_po4_h) * x_h
    rho[5] = params.mu_max_h * f_t_h * f_ph_h * (1.0 - kn_knh_h) * min(m_ss, m_o2_h, m_no3_h, m_po4_h) * x_h
    rho[6] = params.b_max_r_h * f_t_h * f_ph_h * m_o2_h * x_h
    rho[7] = params.mu_max_h * params.eta_anox * f_t_h * f_ph_h * anox * min(m_ss, m_no2_h, m_po4_h) * x_h
    rho[8] = params.mu_max_h * params.eta_anox * f_t_h * f_ph_h * anox * min(m_ss, m_no3_h, m_po4_h) * x_h
    rho[9] = params.b_max_r_h * params.eta_anox * f_t_h * f_ph_h * anox * min(m_no2_h, m_no3_h) * x_h

    # --- Hydrolysis & decay ρ11–ρ13 ---
    if x_h > 0.0:
        ratio = x_s / x_h
        sat_hyd = ratio / (params.k_hyd + ratio)
    else:
        sat_hyd = 0.0
    rho[10] = params.mu_hyd * theta_hyd_T * f_ph_h * sat_hyd * x_h

    rho[11] = params.mu_a * theta_amm_T * f_ph_h * _monod(float(conc[S_ND]), params.k_a) * x_h
    rho[12] = params.b_max_d_h * theta_h_T * f_ph_h * x_h

    # --- AOB ρ14–ρ16, NOB ρ17–ρ19 ---
    m_nh_aob = _monod(float(conc[S_NH]), params.k_n_aob)
    m_o2_aob = _monod(s_o2, params.k_o_aob)
    m_ic_aob = _monod(float(conc[S_IC]), params.k_c_aob)
    m_po4_aob = _monod(float(conc[S_PO4]), params.k_p_aob)
    rho[13] = params.mu_max_aob * f_t_aob * f_ph_aob * min(m_nh_aob, m_o2_aob, m_ic_aob, m_po4_aob) * x_aob
    rho[14] = params.b_max_r_aob * f_t_aob * f_ph_aob * m_o2_aob * x_aob
    rho[15] = params.b_max_d_aob * theta_aob_T * f_ph_aob * x_aob

    m_no2_nob = _monod(float(conc[S_NO2]), params.k_no2_nob)
    m_o2_nob = _monod(s_o2, params.k_o_nob)
    m_ic_nob = _monod(float(conc[S_IC]), params.k_c_nob)
    m_po4_nob = _monod(float(conc[S_PO4]), params.k_p_nob)
    rho[16] = params.mu_max_nob * f_t_nob * f_ph_nob * min(m_no2_nob, m_o2_nob, m_ic_nob, m_po4_nob) * x_nob
    rho[17] = params.b_max_r_nob * f_t_nob * f_ph_nob * m_o2_nob * x_nob
    rho[18] = params.b_max_d_nob * theta_nob_T * f_ph_nob * x_nob

    return rho
