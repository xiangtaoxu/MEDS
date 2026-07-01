#!/usr/bin/env python3
"""Reproduce Slot & Winter (2017) Figure 2 by calling the MEDS leaf-physiology model FROM PYTHON.

The species parameters (the paper's Table 2), the Medlyn-2002 -> model temperature-response
conversion, the VPD(T) relationship and the leaf-temperature sweep all live HERE, in Python, while
the actual photosynthesis kernels are the SAME compiled Fortran the demographic engine uses, called
through post_proc/meds_leaf.py -> libmeds_leaf_c. It is a showcase of MEDS's modularity: the model
lives in Fortran, but no parameters are hard-coded there — the whole experiment is a Python script.

Prerequisites:
    cmake -S . -B build-py -DCMAKE_Fortran_COMPILER=ifx -DMEDS_BUILD_PYLIB=ON -DMEDS_ENABLE_IO=OFF
    cmake --build build-py --target meds_leaf_c
    source /opt/intel/oneapi/setvars.sh        # put the Intel runtime on LD_LIBRARY_PATH

Run (writes CSVs + the figure into this folder):
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
from meds.leaf import (gas_exchange, peaked, arrhenius, make_params,     # noqa: E402
                       Pathway, Stomata, TempResponse, Colimitation)
from plot_slot2017_fig2 import plot, SPECIES                             # noqa: E402

#----- Constants (match meds_constants + the Fortran driver exactly). ----------------------#
R = 8.314462618
T_KELVIN = 273.15
PAR, CA, PRESSURE = 1500.0, 400.0, 101325.0       # paper's A-Ci conditions
EAIR_KPA, G1, G0, RD_FRAC, EA_RD = 2.8, 4.0, 0.01, 0.005, 46390.0
TC = np.arange(25.0, 42.0 + 1e-9, 0.5)            # leaf temperature sweep [degC]

#----- Slot & Winter (2017) Plant Cell Environ. 40:3055, Table 2. Order: F. insipida,        #
#       L. speciosa, C. longifolium, G. madruno.  Each entry: (TOpt[degC], kOpt, Ha, Hd[kJ/mol]).#
#       M4 = four free parameters (solid lines); M3 = Hd fixed at 200 (dashed lines). --------#
VC_M4 = [(36.0, 218, 77, 1049), (39.7, 346, 79, 2975), (32.9, 159, 349, 450), (37.1, 73, 69, 875)]
VC_M3 = [(36.3, 194, 121, 200), (40.3, 278, 107, 200), (33.5, 143, 108, 200), (38.9, 73, 109, 200)]
JM_M4 = [(34.5, 214, 48, 600), (37.5, 226, 65, 610), (33.5, 155, 98, 467), (35.3, 82, 250, 266)]
JM_M3 = [(33.1, 202, 93, 200), (37.0, 205, 98, 200), (31.8, 141, 78, 200), (35.6, 83, 51, 200)]


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
    """Leaf-to-air VPD [Pa] from leaf temperature: saturation deficit above a fixed air e."""
    esat_kpa = 0.61078 * math.exp(17.27 * tc / (tc + 237.3))   # Tetens
    return max((esat_kpa - EAIR_KPA) * 1e3, 100.0)


def run(prefix):
    for isp, sp in enumerate(SPECIES):
        v4 = to_model_form(*VC_M4[isp]); v3 = to_model_form(*VC_M3[isp])
        j4 = to_model_form(*JM_M4[isp]); j3 = to_model_form(*JM_M3[isp])
        k25v4, eav4, hdv4, dsv4 = v4
        k25j4, eaj4, hdj4, dsj4 = j4
        #----- One flat parameter set for the coupled solve (M4 biochemistry). --------------#
        p = make_params(pathway=Pathway.C3, vcmax25=k25v4, jmax25=k25j4,
                        ea_vcmax=eav4, hd_vcmax=hdv4, ds_vcmax=dsv4,
                        ea_jmax=eaj4, hd_jmax=hdj4, ds_jmax=dsj4,
                        rd25=RD_FRAC * k25v4, ea_rd=EA_RD, hd_rd=1e9, ds_rd=490.0,
                        tpu25=1e6, kp25=0.0, g0=G0, g1=G1, d0=1500.0)
        with open(f"{prefix}_{sp}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["tleaf_c", "vcmax_m4", "vcmax_m3", "jmax_m4", "jmax_m3", "rlight", "gs", "anet"])
            for tc in TC:
                tk = tc + T_KELVIN
                vcmax_m4 = peaked(*v4, tk)
                vcmax_m3 = peaked(*v3, tk)
                jmax_m4 = peaked(*j4, tk)
                jmax_m3 = peaked(*j3, tk)
                rlight = arrhenius(RD_FRAC * k25v4, EA_RD, tk)
                flux = gas_exchange(par=PAR, leaf_temp=tk, vpd=vpd_from_temp(tc), ca=CA, params=p,
                                    pressure=PRESSURE, psi_leaf=0.0, gb=0.0,
                                    stomata=Stomata.MEDLYN, temp_response=TempResponse.PEAKED,
                                    colimitation=Colimitation.QUADRATIC, boundary_layer=False)
                w.writerow([f"{tc:.2f}"] + [f"{x:.7e}" for x in
                           (vcmax_m4, vcmax_m3, jmax_m4, jmax_m3, rlight, flux.gs, flux.a_net)])
        print(f"  {sp}: Vcmax25={k25v4:.1f} Jmax25={k25j4:.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-p", "--prefix", default=os.path.join(HERE, "slot2017", "slot2017"),
                    help="output CSV prefix")
    ap.add_argument("-o", "--out", default=os.path.join(HERE, "slot2017_fig2.png"),
                    help="output figure PNG")
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.prefix), exist_ok=True)
    print("Reproducing Slot & Winter (2017) Fig. 2 via the MEDS model (Python -> libmeds_leaf_c):")
    run(args.prefix)
    plot(args.prefix, args.out)


if __name__ == "__main__":
    main()
