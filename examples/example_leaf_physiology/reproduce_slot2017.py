#!/usr/bin/env python3
"""Reproduce Slot & Winter (2017) Figures 1(b) + 2 as ONE figure, calling the MEDS leaf model FROM PYTHON.

Everything lives HERE, in Python (species parameters from the paper's Table 2, the humidity assumption
and the sweeps), while the photosynthesis kernels are the SAME compiled Fortran the demographic engine
uses, reached through the `meds.leaf` package -> libmeds_plant_c. A showcase of MEDS's modularity: the
model lives in Fortran, but no parameters are hard-coded there.

One consolidated figure (slot2017.png):
  * LEFT  -- Fig 1(b)-style A-Ci demand curve for F. insipida, using the paper's CORRECTED in-situ
    Vcmax = 161 / Jmax = 238 DIRECTLY (no capacity temperature-correction). It composes the model
    kernels -- kinetics via arrhenius, J via electron_transport_j, then assimilation_demand_c3 with a SHARP
    minimum -- so the net "limiting rate" A coincides with the lower of the RuBP-carboxylation (Ac) and
    RuBP-regeneration (Aj) net curves, exactly as in the paper.
  * RIGHT -- Fig 2 as five stacked leaf-temperature response panels (Vcmax, Jmax, gs, Anet, Rlight),
    species distinguished by colour, via the coupled solver.

The stomatal/humidity settings (constant RH, Medlyn g1/g0) are tuned so the coupled gs and Anet peak
near ambient temperature as observed in the paper's Fig 2 -- see the README for the exploration.

Prerequisites:
    cmake -S . -B build-py -DCMAKE_Fortran_COMPILER=ifx -DMEDS_BUILD_PYLIB=ON -DMEDS_ENABLE_IO=OFF
    cmake --build build-py --target meds_plant_c
    source /opt/intel/oneapi/setvars.sh        # put the Fortran runtime on LD_LIBRARY_PATH

Run (writes CSVs + slot2017.png into this folder):
    python examples/example_leaf_physiology/reproduce_slot2017.py
"""
import os
import sys
import csv
import math
import argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "python"))     # import meds.leaf from source (no pip install needed)
sys.path.insert(0, os.path.join(ROOT, "post_proc"))
from meds.leaf import (gas_exchange, assimilation_demand_c3, electron_transport_j,     # noqa: E402
                       peaked, arrhenius, make_params,
                       Pathway, Stomata, TempResponse, Colimitation)
from plot_slot2017 import plot_combined, SPECIES                                # noqa: E402

#----- Constants (match meds_constants). ---------------------------------------------------#
R, T_KELVIN = 8.314462618, 273.15
PAR, CA, PRESSURE = 1500.0, 400.0, 101325.0       # paper's A-Ci light + reference CO2
RD_FRAC, EA_RD = 0.005, 46390.0                   # Rd25 = 0.5% of Vcmax25; Rd activation energy
#----- Humidity + stomatal settings tuned to the paper's coordinated gs/Anet temperature optima. -#
#       Using a CONSTANT RELATIVE HUMIDITY (the humid tropical environment) rather than a fixed air  #
#       vapour pressure keeps leaf VPD moderate as temperature rises (0.95 -> 2.5 kPa over 25-42 C   #
#       instead of 0.4 -> 5.4 kPa), so the semi-empirical Medlyn gs tracks A and peaks near ambient   #
#       -- reproducing the paper's central finding. (Exploration + rationale: see the README.) -------#
REL_HUMIDITY = 0.70
G1_MEDLYN, G0_STOM = 4.0, 0.02
TC = np.arange(25.0, 42.0 + 1e-9, 0.5)            # leaf-temperature sweep [degC]
#----- A-Ci panel (Fig 1b): the paper's CORRECTED in-situ capacities for F. insipida, used directly. -#
VCMAX_ACI, JMAX_ACI = 161.0, 238.0                # [umol/m2/s] Slot & Winter (2017) Fig 1b inset

#----- Slot & Winter (2017) Table 2, corrected 4-parameter peaked fits. Order: F. insipida,   #
#       L. speciosa, C. longifolium, G. madruno.  Each entry: (TOpt[degC], kOpt, Ha, Hd[kJ/mol]).#
VC_M4 = [(36.0, 218, 77, 1049), (39.7, 346, 79, 2975), (32.9, 159, 349, 450), (37.1, 73, 69, 875)]
JM_M4 = [(34.5, 214, 48, 600), (37.5, 226, 65, 610), (33.5, 155, 98, 467), (35.3, 82, 250, 266)]


def to_model_form(topt_c, kopt, ha_kj, hd_kj):
    """Convert the paper's peaked form (TOpt, kOpt, Ha, Hd) to the model's (k25, Ea, Hd, dS).

    TOpt = Hd/(dS - R ln(Ha/(Hd-Ha))) is inverted for dS; k25 is set so the model's peaked curve
    passes through kOpt at TOpt (using the model's own peaked function for the shape)."""
    toptk = topt_c + T_KELVIN
    ea, hd = ha_kj * 1e3, hd_kj * 1e3
    ds = hd / toptk + R * math.log(ea / (hd - ea))
    k25 = kopt / peaked(1.0, ea, hd, ds, toptk)          # shape at TOpt via the model
    return k25, ea, hd, ds


def vpd_from_temp(tc):
    """Leaf VPD [Pa] at constant relative humidity: esat(T) * (1 - RH). Tetens saturation."""
    esat_kpa = 0.61078 * math.exp(17.27 * tc / (tc + 237.3))
    return max(esat_kpa * (1.0 - REL_HUMIDITY) * 1e3, 100.0)


def species_params(isp):
    """Corrected-M4 photosynthetic + stomatal Params for Slot species `isp`, plus its (Vc, Jm)
    peaked tuples (k25, Ea, Hd, dS) for drawing the capacity curves."""
    v = to_model_form(*VC_M4[isp])
    j = to_model_form(*JM_M4[isp])
    p = make_params(pathway=Pathway.C3, vcmax25=v[0], jmax25=j[0],
                    ea_vcmax=v[1], hd_vcmax=v[2], ds_vcmax=v[3],
                    ea_jmax=j[1], hd_jmax=j[2], ds_jmax=j[3],
                    rd25=RD_FRAC * v[0], ea_rd=EA_RD, hd_rd=1e9, ds_rd=490.0,
                    tpu25=1e6, kp25=0.0, g0=G0_STOM, g1=G1_MEDLYN, d0=1500.0)
    return p, v, j


def run_temperature(prefix):
    """Fig 2 data: leaf-temperature sweep of Vcmax, Jmax, Rlight and the coupled gs / Anet per species."""
    for isp, sp in enumerate(SPECIES):
        p, v, j = species_params(isp)
        with open(f"{prefix}_{sp}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["tleaf_c", "vcmax", "jmax", "rlight", "gs", "anet"])
            for tc in TC:
                tk = tc + T_KELVIN
                flux = gas_exchange(par=PAR, leaf_temp=tk, vpd=vpd_from_temp(tc), ca=CA, params=p,
                                    pressure=PRESSURE, stomata=Stomata.MEDLYN,
                                    temp_response=TempResponse.PEAKED,
                                    colimitation=Colimitation.QUADRATIC)
                w.writerow([f"{tc:.2f}"] + [f"{x:.7e}" for x in
                           (peaked(*v, tk), peaked(*j, tk), arrhenius(RD_FRAC * v[0], EA_RD, tk),
                            flux.gs, flux.A_net)])
        print(f"  {sp}: Vcmax25={v[0]:.1f} Jmax25={j[0]:.1f}")


def run_aci(prefix, ci_max=1400.0, npts=90):
    """Fig 1(b) data: A-Ci demand curve for F. insipida at 25 C, using Vcmax = 161 / Jmax = 238
    DIRECTLY (no capacity temperature-correction). Composes the model kernels: mole-fraction kinetics
    via arrhenius, J via electron_transport_j, then assimilation_demand_c3 with a SHARP minimum. Writes the
    NET rates (minus Rd) so the limiting rate coincides with the lower of the Ac / Aj net curves. The
    sweep starts near Gamma* (~42 ppm at 25 C), the compensation region."""
    t_leaf = 25.0 + T_KELVIN
    p = make_params(pathway=Pathway.C3, vcmax25=VCMAX_ACI, jmax25=JMAX_ACI)   # Bernacchi kinetics defaults
    kc    = arrhenius(p.kc25,    p.ea_kc,    t_leaf) / PRESSURE * 1e6         # Pa -> mole fraction [ppm]
    ko    = arrhenius(p.ko25,    p.ea_ko,    t_leaf) / PRESSURE * 1e6
    gstar = arrhenius(p.gstar25, p.ea_gstar, t_leaf) / PRESSURE * 1e6
    o2    = p.o2_mol_frac * 1e6
    jrate = electron_transport_j(PAR, JMAX_ACI, absorptance=p.absorptance, phi_psii=p.phi_psii,
                                 theta=p.theta_j)
    rd    = RD_FRAC * VCMAX_ACI                                              # Rd at 25 C
    with open(f"{prefix}_aci.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ci", "ac", "aj", "anet"])                              # NET rates (gross - Rd)
        for ci in np.linspace(40.0, ci_max, npts):
            r = assimilation_demand_c3(ci, VCMAX_ACI, jrate, tpu=p.tpu25, gstar=gstar, kc=kc, ko=ko, o2=o2,
                                colimitation=Colimitation.MINIMUM, theta=p.theta_j)
            w.writerow([f"{ci:.2f}"] + [f"{x:.7e}" for x in (r.Ac - rd, r.Aj - rd, r.A_gross - rd)])
    print(f"  A-Ci (F. insipida): Vcmax={VCMAX_ACI:.0f} Jmax={JMAX_ACI:.0f} J={jrate:.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-p", "--prefix", default=os.path.join(HERE, "slot2017", "slot2017"),
                    help="output CSV prefix")
    ap.add_argument("-o", "--out", default=os.path.join(HERE, "slot2017.png"),
                    help="output figure PNG")
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.prefix), exist_ok=True)
    print("Reproducing Slot & Winter (2017) Figs 1b + 2 via the MEDS model (Python -> libmeds_plant_c):")
    run_aci(args.prefix)
    run_temperature(args.prefix)
    plot_combined(args.prefix, args.out)


if __name__ == "__main__":
    main()
