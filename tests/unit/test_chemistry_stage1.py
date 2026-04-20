"""Stage 1 tests: molar conversions, reference K_a, van't Hoff (SI.6.4)."""

from __future__ import annotations

import math

from bioprocess_twin.models.chemistry import (
    R_GAS,
    T_REF_K,
    default_dissociation_constants_ref_molar,
    default_dissociation_enthalpy_j_per_mol,
    ka_at_T,
    kw_at_T,
    mol_c_per_m3_carbon_from_s_ic,
    mol_n_per_m3_nitrogen_from_s_nh,
    mol_n_per_m3_nitrogen_from_s_no2,
    mol_n_per_m3_nitrogen_from_s_no3,
    mol_p_per_m3_phosphorus_from_s_po4,
    scale_dissociation_constants_at_t,
)


def test_ka_at_t_identity_at_t_ref() -> None:
    ka_ref = 1.0e-9
    dh = 50_000.0
    t_celsius = T_REF_K - 273.15
    assert math.isclose(ka_at_T(ka_ref, dh, t_celsius), ka_ref, rel_tol=0.0, abs_tol=1e-15)


def test_ka_at_t_hand_check_warmer_temperature() -> None:
    ka_ref = 1.0e-7
    delta_h = 50_000.0
    t_cold = 15.0
    t_hot = 35.0
    k_cold = ka_at_T(ka_ref, delta_h, t_cold)
    k_hot = ka_at_T(ka_ref, delta_h, t_hot)
    t_k_cold = t_cold + 273.15
    t_k_hot = t_hot + 273.15
    expected_ratio = math.exp(delta_h / R_GAS * (1.0 / t_k_cold - 1.0 / t_k_hot))
    assert math.isclose(k_hot / k_cold, expected_ratio, rel_tol=1e-12)
    assert k_hot > k_cold


def test_default_nh4_ka_matches_si6_table() -> None:
    ref = default_dissociation_constants_ref_molar()
    assert math.isclose(ref.ka_nh4, 5.62e-10, rel_tol=0.0, abs_tol=1e-15)


def test_default_hno3_ka_matches_si6_table() -> None:
    ref = default_dissociation_constants_ref_molar()
    assert math.isclose(ref.ka_hno3, 4.37e1, rel_tol=0.0, abs_tol=1e-15)


def test_default_co2_ka_matches_si6_table() -> None:
    ref = default_dissociation_constants_ref_molar()
    assert math.isclose(ref.ka_co2, 4.27e-7, rel_tol=0.0, abs_tol=1e-15)


def test_mol_c_from_s_ic() -> None:
    assert math.isclose(mol_c_per_m3_carbon_from_s_ic(24.0), 2.0, rel_tol=0.0, abs_tol=1e-15)


def test_mol_n_from_s_nh() -> None:
    assert math.isclose(mol_n_per_m3_nitrogen_from_s_nh(14.0), 1.0, rel_tol=0.0, abs_tol=1e-15)


def test_mol_n_from_s_no2_and_s_no3() -> None:
    assert math.isclose(mol_n_per_m3_nitrogen_from_s_no2(7.0), 0.5, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(mol_n_per_m3_nitrogen_from_s_no3(28.0), 2.0, rel_tol=0.0, abs_tol=1e-15)


def test_mol_p_from_s_po4() -> None:
    assert math.isclose(mol_p_per_m3_phosphorus_from_s_po4(62.0), 2.0, rel_tol=0.0, abs_tol=1e-15)


def test_scale_bundle_matches_individual_ka_at_t() -> None:
    ref = default_dissociation_constants_ref_molar()
    dh = default_dissociation_enthalpy_j_per_mol()
    t_c = 20.0
    scaled = scale_dissociation_constants_at_t(ref, dh, t_c)
    assert math.isclose(scaled.ka_co2, ka_at_T(ref.ka_co2, dh.dh_co2, t_c), rel_tol=1e-15)
    assert math.isclose(scaled.kw, ka_at_T(ref.kw, dh.dh_kw, t_c), rel_tol=1e-15)


def test_zero_delta_h_leaves_ka_unchanged_at_any_t() -> None:
    ka_ref = 2.0e-5
    assert math.isclose(ka_at_T(ka_ref, 0.0, 10.0), ka_ref, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(ka_at_T(ka_ref, 0.0, 40.0), ka_ref, rel_tol=0.0, abs_tol=1e-15)


def test_kw_at_t_via_scale_matches_ka_at_t() -> None:
    ref = default_dissociation_constants_ref_molar()
    dh = default_dissociation_enthalpy_j_per_mol()
    t_c = 22.0
    assert math.isclose(kw_at_T(ref.kw, dh.dh_kw, t_c), ka_at_T(ref.kw, dh.dh_kw, t_c), rel_tol=0.0, abs_tol=1e-15)


def test_hno3_enthalpy_zero_scales_to_identity() -> None:
    ref = default_dissociation_constants_ref_molar()
    dh = default_dissociation_enthalpy_j_per_mol()
    assert dh.dh_hno3 == 0.0
    scaled = scale_dissociation_constants_at_t(ref, dh, 18.0)
    assert math.isclose(scaled.ka_hno3, ref.ka_hno3, rel_tol=0.0, abs_tol=1e-15)
