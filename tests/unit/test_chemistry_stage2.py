"""Stage 2 tests: SI.6.1 speciation (mass closure, limits, hand carbonate, T with dh=0)."""

from __future__ import annotations

import math

from bioprocess_twin.models.chemistry import (
    AlbaDissociationEnthalpy,
    SpeciationTotals,
    default_dissociation_constants_ref_molar,
    mol_c_per_m3_carbon_from_s_ic,
    mol_n_per_m3_nitrogen_from_s_nh,
    mol_n_per_m3_nitrogen_from_s_no2,
    mol_n_per_m3_nitrogen_from_s_no3,
    mol_p_per_m3_phosphorus_from_s_po4,
    scale_dissociation_constants_at_t,
    speciate_ammonia,
    speciate_aqueous,
    speciate_carbonate,
    speciate_from_alba_totals,
)


def test_mass_closure_all_subsystems() -> None:
    """speciate_aqueous preserves each subsystem total (NH, NO2, NO3, IC, P pools sum to c_tot)."""
    k = default_dissociation_constants_ref_molar()
    h = 1.0e-6
    totals = SpeciationTotals(
        c_tot_nh3=0.01,
        c_tot_no2=0.002,
        c_tot_no3=0.003,
        c_tot_ic=0.05,
        c_tot_po4=0.001,
    )
    sp = speciate_aqueous(h, totals, k)
    assert math.isclose(sp.nh3 + sp.nh4, totals.c_tot_nh3, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(sp.no2 + sp.hno2, totals.c_tot_no2, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(sp.no3 + sp.hno3, totals.c_tot_no3, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(sp.co2 + sp.hco3 + sp.co3, totals.c_tot_ic, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        sp.h3po4 + sp.h2po4 + sp.hpo4 + sp.po4,
        totals.c_tot_po4,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_ammonia_low_ph_mostly_nh4() -> None:
    """High [H+]: almost all ammonia nitrogen is NH4+ (acidic form dominates)."""
    k = default_dissociation_constants_ref_molar()
    c_tot = 1.0
    h = 1.0e-2
    nh3, nh4 = speciate_ammonia(h, k.ka_nh4, c_tot)
    assert nh4 > 0.99 * c_tot
    assert nh3 < 0.02 * c_tot


def test_ammonia_high_ph_mostly_nh3() -> None:
    """Very low [H+]: almost all ammonia nitrogen is NH3 (base form dominates)."""
    k = default_dissociation_constants_ref_molar()
    c_tot = 1.0
    h = 1.0e-12
    nh3, nh4 = speciate_ammonia(h, k.ka_nh4, c_tot)
    assert nh3 > 0.95 * c_tot
    assert nh4 < 0.05 * c_tot


def test_carbonate_low_ph_mostly_co2() -> None:
    """High [H+]: inorganic carbon is mostly CO2(aq); carbonate ion fraction stays small."""
    k = default_dissociation_constants_ref_molar()
    c_tot = 0.02
    h = 1.0e-2
    co2, hco3, co3 = speciate_carbonate(h, k.ka_co2, k.ka_hco3, c_tot)
    assert co2 > 0.85 * c_tot
    assert co3 < 0.05 * c_tot


def test_carbonate_high_ph_mostly_co3() -> None:
    """Very low [H+]: CO3^2- is a large share of total IC; CO2 fraction is limited."""
    k = default_dissociation_constants_ref_molar()
    c_tot = 0.02
    h = 1.0e-11
    co2, hco3, co3 = speciate_carbonate(h, k.ka_co2, k.ka_hco3, c_tot)
    assert co3 > 0.5 * c_tot
    assert co2 < 0.3 * c_tot


def test_carbonate_hand_calculation() -> None:
    """CO2, HCO3-, CO3^2- from speciate_carbonate match SI.6.1 row 8–9 formulas (hand recomputed)."""
    k = default_dissociation_constants_ref_molar()
    c_tot = 1.0
    h = 1.0e-5
    co2, hco3, co3 = speciate_carbonate(h, k.ka_co2, k.ka_hco3, c_tot)
    denom_co2 = 1.0 + k.ka_co2 * 1.0e3 / h + k.ka_co2 * k.ka_hco3 * 1.0e6 / (h * h)
    expected_co2 = c_tot / denom_co2
    denom_hco3 = 1.0 + h / (k.ka_hco3 * 1.0e3) + k.ka_hco3 * 1.0e3 / h
    expected_hco3 = c_tot / denom_hco3
    expected_co3 = c_tot - expected_co2 - expected_hco3
    assert math.isclose(co2, expected_co2, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(hco3, expected_hco3, rel_tol=0.0, abs_tol=1e-15)
    assert math.isclose(co3, expected_co3, rel_tol=0.0, abs_tol=1e-15)


def test_speciate_from_alba_totals_matches_speciate_aqueous() -> None:
    """Wrapper speciate_from_alba_totals agrees with speciate_aqueous on the same mol m^-3 totals."""
    k = default_dissociation_constants_ref_molar()
    h = 3.0e-7
    s_ic, s_nh, s_no2, s_no3, s_po4 = 24.0, 1.4, 0.0, 2.8, 3.1
    totals = SpeciationTotals(
        c_tot_nh3=mol_n_per_m3_nitrogen_from_s_nh(s_nh),
        c_tot_no2=mol_n_per_m3_nitrogen_from_s_no2(s_no2),
        c_tot_no3=mol_n_per_m3_nitrogen_from_s_no3(s_no3),
        c_tot_ic=mol_c_per_m3_carbon_from_s_ic(s_ic),
        c_tot_po4=mol_p_per_m3_phosphorus_from_s_po4(s_po4),
    )
    a = speciate_aqueous(h, totals, k)
    b = speciate_from_alba_totals(h, s_ic, s_nh, s_no2, s_no3, s_po4, k)
    assert a == b


def test_zero_delta_h_temperature_does_not_change_speciation() -> None:
    """With all ΔH = 0, scaled K at 12 °C vs 37 °C yields identical speciation."""
    ref = default_dissociation_constants_ref_molar()
    dh = AlbaDissociationEnthalpy(
        dh_co2=0.0,
        dh_hco3=0.0,
        dh_nh4=0.0,
        dh_hno2=0.0,
        dh_hno3=0.0,
        dh_h3po4=0.0,
        dh_h2po4=0.0,
        dh_hpo4=0.0,
        dh_kw=0.0,
    )
    totals = SpeciationTotals(
        c_tot_nh3=0.5e-2,
        c_tot_no2=1.0e-4,
        c_tot_no3=2.0e-3,
        c_tot_ic=0.04,
        c_tot_po4=5.0e-4,
    )
    h = 2.0e-7
    k_cold = scale_dissociation_constants_at_t(ref, dh, 12.0)
    k_hot = scale_dissociation_constants_at_t(ref, dh, 37.0)
    sp_c = speciate_aqueous(h, totals, k_cold)
    sp_h = speciate_aqueous(h, totals, k_hot)
    assert sp_c == sp_h
