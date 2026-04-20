"""
ALBA aqueous chemistry: reference dissociation constants and van't Hoff scaling.

Reference K_a and K_w in mol/L are taken from Table SI.6.1 in
docs/supporting_informations/SI.6 (third column, header 293 K).

Equation SI.6.4 in the same file uses the van't Hoff form with T in Celsius and T_ref in the
text given as 298.15 K; the tabulated K values are at 293 K. This module uses T_REF_K = 293.15
so K_ref in ka_at_T matches the table temperature (same equation shape as SI.6.4).

default_dissociation_enthalpy_j_per_mol uses evaluated standard reaction enthalpies (infinite
dilution, ~298.15 K, 1 bar) form various resources. Mixing these ΔH° with K_ref anchored at 
T_REF_K = 293.15 K is a small inconsistency (~5 K); for full alignment either move K_ref to 298.15 K or 
re-fit ΔH° to the SI.6.1 table temperature.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# SI.6.4: gas constant [J K^-1 mol^-1]
R_GAS = 8.314462618

# Temperature of K_A in Table SI.6.1 (293 K) [K]; used as T_ref in ka_at_T with tabulated K
T_REF_K = 293.15

# Molar masses for ALBA state totals [g mol^-1] (SI.6 uses 12, 14, 31)
M_C = 12.0
M_N = 14.0
M_P = 31.0


@dataclass(frozen=True, slots=True)
class AlbaDissociationConstantsRef:
    """Dissociation constants at T_REF_K (298.15 K), in mol/L. Field names match SI.6.1."""

    ka_co2: float
    ka_hco3: float
    ka_nh4: float
    ka_hno2: float
    ka_hno3: float
    ka_h3po4: float
    ka_h2po4: float
    ka_hpo4: float
    kw: float


@dataclass(frozen=True, slots=True)
class AlbaDissociationEnthalpy:
    """
    Standard reaction enthalpy Delta H in J/mol per equilibrium for SI.6.4.

    Same keys as AlbaDissociationConstantsRef; dh_kw is for water ionization.
    """

    dh_co2: float
    dh_hco3: float
    dh_nh4: float
    dh_hno2: float
    dh_hno3: float
    dh_h3po4: float
    dh_h2po4: float
    dh_hpo4: float
    dh_kw: float


def default_dissociation_constants_ref_molar() -> AlbaDissociationConstantsRef:
    """
    Build K_a,ref from docs/MATH_MODEL.md section 1.2.7 (pK_a at 298.15 K).

    Returns all values in mol/L. kw is 1e-14 at 25 C (common convention).
    """
    return AlbaDissociationConstantsRef(
        ka_co2=10.0**-6.37,
        ka_hco3=10.0**-10.33,
        ka_nh4=10.0**-9.25,
        ka_hno2=10.0**-3.35,
        ka_hno3=10.0**1.64,
        ka_h3po4=10.0**-2.14,
        ka_h2po4=10.0**-7.21,
        ka_hpo4=10.0**-12.67,
        kw=1.0e-14,
    )


def default_dissociation_enthalpy_j_per_mol() -> AlbaDissociationEnthalpy:
    """
    Standard reaction enthalpy ΔH° in J/mol for van't Hoff (SI.6.4), K_a = [products]/[acid].

    Values follow the consolidated recommendations in the aqueous-equilibria thermodynamics
    note (CODATA / ATcT / NEA-TDB–based), infinite dilution, 298.15 K, 1 bar. dh_co2 is for
    the lumped CO2(aq)/H2CO3* first dissociation as in bioprocess models; dh_hno3 is 0 under
    the strong-acid convention (undissociated HNO3(aq) ≡ ions at infinite dilution).
    """
    return AlbaDissociationEnthalpy(
        dh_co2=9_155.0,
        dh_hco3=14_700.0,
        dh_nh4=52_201.0,
        dh_hno2=11_400.0,
        dh_hno3=0.0,
        dh_h3po4=-8_480.0,
        dh_h2po4=3_600.0,
        dh_hpo4=14_600.0,
        dh_kw=55_830.0,
    )


def ka_at_T(ka_ref_molar: float, delta_h_j_per_mol: float, t_celsius: float) -> float:
    """
    Scale equilibrium constant K from T_REF_K to t_celsius using SI.6.4.

    ka_ref_molar: K at T_REF_K in mol/L (for K_w use the same formula; units stay consistent with the reference).
    delta_h_j_per_mol: Delta H in J/mol for that reaction.
    t_celsius: temperature in degrees Celsius.

    Returns K at t_celsius in the same concentration units as ka_ref_molar.
    """
    t_k = t_celsius + 273.15
    return ka_ref_molar * math.exp(
        (delta_h_j_per_mol / R_GAS) * (1.0 / T_REF_K - 1.0 / t_k),
    )


def kw_at_T(kw_ref_molar: float, delta_h_kw_j_per_mol: float, t_celsius: float) -> float:
    """Water ionization constant at t_celsius; same math as ka_at_T. kw_ref_molar is often 1e-14 (mol/L)^2 at 25 C."""
    return ka_at_T(kw_ref_molar, delta_h_kw_j_per_mol, t_celsius)


def mol_c_per_m3_carbon_from_s_ic(s_ic_g_per_m3: float) -> float:
    """Moles of inorganic carbon per m3 from S_IC in g C/m3 (SI.6 row 7)."""
    return s_ic_g_per_m3 / M_C


def mol_n_per_m3_nitrogen_from_s_nh(s_nh_g_per_m3: float) -> float:
    """Moles of ammoniacal nitrogen per m3 from S_NH in g N/m3."""
    return s_nh_g_per_m3 / M_N


def mol_n_per_m3_nitrogen_from_s_no2(s_no2_g_per_m3: float) -> float:
    """Moles of nitrite nitrogen per m3 from S_NO2 in g N/m3."""
    return s_no2_g_per_m3 / M_N


def mol_n_per_m3_nitrogen_from_s_no3(s_no3_g_per_m3: float) -> float:
    """Moles of nitrate nitrogen per m3 from S_NO3 in g N/m3."""
    return s_no3_g_per_m3 / M_N


def mol_p_per_m3_phosphorus_from_s_po4(s_po4_g_per_m3: float) -> float:
    """Moles of inorganic phosphorus per m3 from S_PO4 in g P/m3."""
    return s_po4_g_per_m3 / M_P


def scale_dissociation_constants_at_t(
    ref: AlbaDissociationConstantsRef,
    dh: AlbaDissociationEnthalpy,
    t_celsius: float,
) -> AlbaDissociationConstantsRef:
    """Apply ka_at_T / kw_at_T to each field of ref using dh."""
    return AlbaDissociationConstantsRef(
        ka_co2=ka_at_T(ref.ka_co2, dh.dh_co2, t_celsius),
        ka_hco3=ka_at_T(ref.ka_hco3, dh.dh_hco3, t_celsius),
        ka_nh4=ka_at_T(ref.ka_nh4, dh.dh_nh4, t_celsius),
        ka_hno2=ka_at_T(ref.ka_hno2, dh.dh_hno2, t_celsius),
        ka_hno3=ka_at_T(ref.ka_hno3, dh.dh_hno3, t_celsius),
        ka_h3po4=ka_at_T(ref.ka_h3po4, dh.dh_h3po4, t_celsius),
        ka_h2po4=ka_at_T(ref.ka_h2po4, dh.dh_h2po4, t_celsius),
        ka_hpo4=ka_at_T(ref.ka_hpo4, dh.dh_hpo4, t_celsius),
        kw=kw_at_T(ref.kw, dh.dh_kw, t_celsius),
    )


# --- Stage 2: speciation (SI.6.1, rows 2–14); [H+] and all species in mol m^-3; K_a, K_w in mol/L ---


@dataclass(frozen=True, slots=True)
class SpeciationTotals:
    """Total molar concentrations per SI.6.1 mass-balance rows (mol m^-3)."""

    c_tot_nh3: float
    """Total ammonia nitrogen NH3 + NH4+ (mol N m^-3), SI row 1 (as S_NH3/14)."""

    c_tot_no2: float
    """Total nitrite nitrogen NO2- + HNO2 (mol N m^-3), SI row 3."""

    c_tot_no3: float
    """Total nitrate nitrogen NO3- + HNO3 (mol N m^-3), SI row 5."""

    c_tot_ic: float
    """Total inorganic carbon CO2 + HCO3- + CO3^2- (mol C m^-3), SI row 7."""

    c_tot_po4: float
    """Total orthophosphate (all protonation states) (mol P m^-3), SI row 10."""


@dataclass(frozen=True, slots=True)
class AqueousSpeciesMolar:
    """Aqueous species concentrations (mol m^-3) from SI.6.1 rows 2, 4, 6, 8–9, 11–13, 14."""

    nh3: float
    nh4: float
    no2: float
    hno2: float
    no3: float
    hno3: float
    co2: float
    hco3: float
    co3: float
    h3po4: float
    h2po4: float
    hpo4: float
    po4: float
    oh: float


def speciate_ammonia(h_plus_mol_per_m3: float, ka_nh4_molar: float, c_tot_nh_mol_per_m3: float) -> tuple[float, float]:
    """
    SI.6.1 row 2; K_a in mol/L. Returns (NH3, NH4+) in mol m^-3.

    The SI table labels the quotient as NH3, but the same algebraic form is [NH4+] =
    C_T * h / (h + K_a * 10^3) in this unit convention; NH3 is the complement for mass balance.
    """
    denom = 1.0 + ka_nh4_molar * 1.0e3 / h_plus_mol_per_m3
    nh4 = c_tot_nh_mol_per_m3 / denom
    nh3 = c_tot_nh_mol_per_m3 - nh4
    return nh3, nh4


def speciate_nitrite(
    h_plus_mol_per_m3: float,
    ka_hno2_molar: float,
    c_tot_no2_mol_per_m3: float,
) -> tuple[float, float]:
    """SI.6.1 row 4. Returns (NO2-, HNO2) in mol m^-3."""
    denom = 1.0 + ka_hno2_molar * 1.0e3 / h_plus_mol_per_m3
    hno2 = c_tot_no2_mol_per_m3 / denom
    no2 = c_tot_no2_mol_per_m3 - hno2
    return no2, hno2


def speciate_nitrate(
    h_plus_mol_per_m3: float,
    ka_hno3_molar: float,
    c_tot_no3_mol_per_m3: float,
) -> tuple[float, float]:
    """SI.6.1 row 6. Returns (NO3-, HNO3) in mol m^-3."""
    denom = 1.0 + ka_hno3_molar * 1.0e3 / h_plus_mol_per_m3
    hno3 = c_tot_no3_mol_per_m3 / denom
    no3 = c_tot_no3_mol_per_m3 - hno3
    return no3, hno3


def speciate_carbonate(
    h_plus_mol_per_m3: float,
    ka_co2_molar: float,
    ka_hco3_molar: float,
    c_tot_ic_mol_per_m3: float,
) -> tuple[float, float, float]:
    """
    SI.6.1 rows 8–9; CO3^2- from row 7 closure. Returns (CO2, HCO3-, CO3^2-) in mol m^-3.
    """
    h = h_plus_mol_per_m3
    denom_co2 = 1.0 + ka_co2_molar * 1.0e3 / h + ka_co2_molar * ka_hco3_molar * 1.0e6 / (h * h)
    co2 = c_tot_ic_mol_per_m3 / denom_co2
    denom_hco3 = 1.0 + h / (ka_hco3_molar * 1.0e3) + ka_hco3_molar * 1.0e3 / h
    hco3 = c_tot_ic_mol_per_m3 / denom_hco3
    co3 = c_tot_ic_mol_per_m3 - co2 - hco3
    return co2, hco3, co3


def speciate_phosphate(
    h_plus_mol_per_m3: float,
    ka_h3po4_molar: float,
    ka_h2po4_molar: float,
    ka_hpo4_molar: float,
    c_tot_p_mol_per_m3: float,
) -> tuple[float, float, float, float]:
    """
    SI.6.1 rows 11–13; PO4^3- from row 10 closure.
    Returns (H3PO4, H2PO4-, HPO4^2-, PO4^3-) in mol m^-3.
    """
    h = h_plus_mol_per_m3
    denom_h3po4 = (
        1.0
        + ka_h3po4_molar * 1.0e3 / h
        + ka_h3po4_molar * ka_h2po4_molar * 1.0e6 / (h * h)
        + ka_h3po4_molar * ka_h2po4_molar * ka_hpo4_molar * 1.0e9 / (h * h * h)
    )
    h3po4 = c_tot_p_mol_per_m3 / denom_h3po4
    denom_h2po4 = (
        1.0
        + h / (ka_h2po4_molar * 1.0e3)
        + ka_h2po4_molar * 1.0e3 / h
        + ka_h2po4_molar * ka_hpo4_molar * 1.0e6 / (h * h)
    )
    h2po4 = c_tot_p_mol_per_m3 / denom_h2po4
    denom_hpo4 = (
        1.0
        + (h * h) / (ka_h2po4_molar * ka_hpo4_molar * 1.0e6)
        + h / (ka_hpo4_molar * 1.0e3)
        + ka_hpo4_molar * 1.0e3 / h
    )
    hpo4 = c_tot_p_mol_per_m3 / denom_hpo4
    po4 = c_tot_p_mol_per_m3 - h3po4 - h2po4 - hpo4
    return h3po4, h2po4, hpo4, po4


def speciate_water(h_plus_mol_per_m3: float, kw_molar: float) -> float:
    """SI.6.1 row 14. Returns OH- in mol m^-3."""
    return kw_molar * 1.0e3 / h_plus_mol_per_m3


def speciate_aqueous(
    h_plus_mol_per_m3: float,
    totals: SpeciationTotals,
    k: AlbaDissociationConstantsRef,
) -> AqueousSpeciesMolar:
    """
    Full aqueous speciation for SI.6.1 rows 2–14 (given [H+] and totals; no charge balance).
    """
    nh3, nh4 = speciate_ammonia(h_plus_mol_per_m3, k.ka_nh4, totals.c_tot_nh3)
    no2, hno2 = speciate_nitrite(h_plus_mol_per_m3, k.ka_hno2, totals.c_tot_no2)
    no3, hno3 = speciate_nitrate(h_plus_mol_per_m3, k.ka_hno3, totals.c_tot_no3)
    co2, hco3, co3 = speciate_carbonate(
        h_plus_mol_per_m3,
        k.ka_co2,
        k.ka_hco3,
        totals.c_tot_ic,
    )
    h3po4, h2po4, hpo4, po4 = speciate_phosphate(
        h_plus_mol_per_m3,
        k.ka_h3po4,
        k.ka_h2po4,
        k.ka_hpo4,
        totals.c_tot_po4,
    )
    oh = speciate_water(h_plus_mol_per_m3, k.kw)
    return AqueousSpeciesMolar(
        nh3=nh3,
        nh4=nh4,
        no2=no2,
        hno2=hno2,
        no3=no3,
        hno3=hno3,
        co2=co2,
        hco3=hco3,
        co3=co3,
        h3po4=h3po4,
        h2po4=h2po4,
        hpo4=hpo4,
        po4=po4,
        oh=oh,
    )


def speciate_from_alba_totals(
    h_plus_mol_per_m3: float,
    s_ic_g_per_m3: float,
    s_nh_g_per_m3: float,
    s_no2_g_per_m3: float,
    s_no3_g_per_m3: float,
    s_po4_g_per_m3: float,
    k: AlbaDissociationConstantsRef,
) -> AqueousSpeciesMolar:
    """
    Speciation from ALBA mass totals (g m^-3) via SI.6.1 /12, /14, /31 conversions.
    """
    totals = SpeciationTotals(
        c_tot_nh3=mol_n_per_m3_nitrogen_from_s_nh(s_nh_g_per_m3),
        c_tot_no2=mol_n_per_m3_nitrogen_from_s_no2(s_no2_g_per_m3),
        c_tot_no3=mol_n_per_m3_nitrogen_from_s_no3(s_no3_g_per_m3),
        c_tot_ic=mol_c_per_m3_carbon_from_s_ic(s_ic_g_per_m3),
        c_tot_po4=mol_p_per_m3_phosphorus_from_s_po4(s_po4_g_per_m3),
    )
    return speciate_aqueous(h_plus_mol_per_m3, totals, k)
