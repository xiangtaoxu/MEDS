#!/usr/bin/env python3
"""Plot leaf gas-exchange response curves written by the meds_leaf_demo driver.

meds_leaf_demo writes four CSVs (<prefix>_aci.csv, _apar.csv, _atemp.csv, _gsvpd.csv). This
script renders them as a 2x2 panel of canonical leaf-physiology diagnostics:

  (A-Ci)   net assimilation vs intercellular CO2, with the Ac/Aj/Ap limitation envelope;
  (A-PAR)  the light response of the fully coupled solve;
  (A-Tleaf) the temperature response (the peaked-Arrhenius thermal optimum);
  (gs-VPD) stomatal closure with rising VPD for the three stomatal models.

Usage:
    python plot_leaf_response.py PREFIX [-o OUT.png]

Requires the `meds` conda environment (numpy, matplotlib). Inputs are plain CSV (no netCDF).
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")            # headless: write a PNG, no display needed
import matplotlib.pyplot as plt


def load(prefix, suffix):
    """Read <prefix>_<suffix>.csv as a structured array with column names from the header.

    Leading '#'-comment lines (e.g. the pathway marker the demo writes) are skipped so the
    next line is taken as the column header."""
    path = f"{prefix}_{suffix}.csv"
    nhdr = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                nhdr += 1
            else:
                break
    return np.genfromtxt(path, delimiter=",", names=True, skip_header=nhdr)


def read_pathway(prefix):
    """Return 'c3' or 'c4' from the '# pathway = ...' comment the demo writes to the A-Ci CSV."""
    try:
        with open(f"{prefix}_aci.csv") as fh:
            for line in fh:
                if line.startswith("#") and "pathway" in line:
                    return line.split("=", 1)[1].strip().lower()
                if not line.startswith("#"):
                    break
    except OSError:
        pass
    return "c3"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="output prefix passed to meds_leaf_demo (e.g. leafdemo)")
    ap.add_argument("-o", "--out", default=None, help="output PNG (default: <prefix>_leaf.png)")
    args = ap.parse_args()
    out = args.out or (args.prefix + "_leaf.png")

    aci = load(args.prefix, "aci")
    apar = load(args.prefix, "apar")
    atemp = load(args.prefix, "atemp")
    gsvpd = load(args.prefix, "gsvpd")
    pathway = read_pathway(args.prefix)

    #----- The third limitation rate is pathway-specific: C3 = TPU/product, C4 = PEPcase CO2. -#
    if pathway == "c4":
        lbl_ac, lbl_aj, lbl_ap = "Ac (Vcmax/Rubisco)", "Aj (light)", "Ap (PEP / CO$_2$)"
    else:
        lbl_ac, lbl_aj, lbl_ap = "Ac (Rubisco)", "Aj (RuBP/light)", "Ap (TPU / product)"

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))

    #----- A-Ci with the Ac/Aj/Ap limitation envelope. ------------------------------------#
    a = ax[0, 0]
    a.plot(aci["ci"], aci["ac"], color="tab:red", lw=1, ls="--", label=lbl_ac)
    a.plot(aci["ci"], aci["aj"], color="tab:blue", lw=1, ls="--", label=lbl_aj)
    a.plot(aci["ci"], aci["ap"], color="tab:green", lw=1, ls="--", label=lbl_ap)
    a.plot(aci["ci"], aci["a_net"], color="black", lw=2, label="A_net")
    a.set_xlabel("Ci [umol/mol]"); a.set_ylabel("A [umol m$^{-2}$ s$^{-1}$]")
    a.set_title("A-Ci demand curve"); a.legend(fontsize=8); a.grid(alpha=0.3)
    a.set_ylim(bottom=min(0.0, float(np.nanmin(aci["a_net"]))))

    #----- A-PAR light response. ----------------------------------------------------------#
    a = ax[0, 1]
    a.plot(apar["par"], apar["a_net"], color="black", lw=2)
    a.set_xlabel("PAR [umol photon m$^{-2}$ s$^{-1}$]"); a.set_ylabel("A_net [umol m$^{-2}$ s$^{-1}$]")
    a.set_title("Light response (coupled solve)"); a.grid(alpha=0.3)

    #----- A-leaf-temperature response (thermal optimum). ---------------------------------#
    a = ax[1, 0]
    a.plot(atemp["tleaf_c"], atemp["a_net"], color="black", lw=2)
    a.set_xlabel("leaf temperature [degC]"); a.set_ylabel("A_net [umol m$^{-2}$ s$^{-1}$]")
    a.set_title("Temperature response"); a.grid(alpha=0.3)

    #----- gs(VPD) for the three stomatal models. -----------------------------------------#
    a = ax[1, 1]
    a.plot(gsvpd["vpd_pa"], gsvpd["gs_leuning"], color="tab:orange", lw=2, label="Leuning")
    a.plot(gsvpd["vpd_pa"], gsvpd["gs_medlyn"], color="tab:purple", lw=2, label="Medlyn")
    a.plot(gsvpd["vpd_pa"], gsvpd["gs_katul"], color="tab:cyan", lw=2, label="Katul")
    a.set_xlabel("VPD [Pa]"); a.set_ylabel("gs [mol H$_2$O m$^{-2}$ s$^{-1}$]")
    a.set_title("Stomatal closure with VPD"); a.legend(fontsize=8); a.grid(alpha=0.3)

    fig.suptitle(f"MEDS leaf gas exchange ({args.prefix})", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
