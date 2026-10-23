#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""check_integrator_boundary.py -- the integrator steps one patch and sees nothing past it.

src/fast_dynamics/integrator advances one patch over one dt_fast from the records the driver hands
it: the configuration, the patch's site, and its static state, forcing, dynamic state, output and
workspace. It never sees the polygon, the cohort and patch blocks, the forcing reader, the output
layer, the calendar, or the driver that walks the patches: those are src/fast_dynamics/driver's.
As with the kernels' rule (check_kernel_state.py), the shared module directory cannot hold this,
so this script checks it on every ctest run. Code in the integrator:
  * never names polygon_t, cohort_block, patch_block or site_block;
  * never uses meds_polygon_state_vars, meds_cohort_state_vars, a met or output module
    (meds_met*, meds_output*), meds_time, or a driver module;
  * takes from meds_site_state_vars only site_surface_t, and from meds_patch_state_vars only the
    per-patch store types.

Usage: check_integrator_boundary.py SRC_DIR   (exit 1 and a list when the integrator reaches past its patch)
"""
import re
import sys
from pathlib import Path

from check_kernel_state import ALLOWED as PATCH_STORE_TYPES, USE, statements

INTEGRATOR_DIR = "fast_dynamics/integrator"
FORBIDDEN_TYPES = ("polygon_t", "cohort_block", "patch_block", "site_block")
FORBIDDEN_MODULES = ("meds_polygon_state_vars", "meds_cohort_state_vars", "meds_time",
                     "meds_fast_dynamics", "meds_fast_state", "meds_fast_forcing", "meds_fast_diag",
                     "meds_polygon", "meds_region", "meds_stepper")
FORBIDDEN_PREFIXES = ("meds_met", "meds_output")
ONLY = {"meds_site_state_vars": {"site_surface_t"}, "meds_patch_state_vars": PATCH_STORE_TYPES}


def violations(src):
    found = []
    for f in sorted((src / INTEGRATOR_DIR).glob("*.f90")):
        for no, stmt in statements(f):
            for t in FORBIDDEN_TYPES:
                if re.search(rf"\b{t}\b", stmt, re.I):
                    found.append(f"{f}:{no}: names {t}: {stmt[:100]}")
            m = USE.match(stmt)
            if not m:
                continue
            mod = m.group(1).lower()
            if mod in FORBIDDEN_MODULES or mod.startswith(FORBIDDEN_PREFIXES):
                found.append(f"{f}:{no}: uses {mod}")
            elif mod in ONLY:
                if not m.group(2):
                    found.append(f"{f}:{no}: uses {mod} without an only-list")
                    continue
                names = [n.strip().split("=>")[-1].strip().lower() for n in m.group(3).split(",") if n.strip()]
                for n in names:
                    if n not in ONLY[mod]:
                        found.append(f"{f}:{no}: takes {n} from {mod}")
    return found


if __name__ == "__main__":
    src = Path(sys.argv[1] if len(sys.argv) > 1 else "src")
    bad = violations(src)
    for line in bad:
        print(line)
    if bad:
        print(f"\n{len(bad)} integrator reference(s) past the patch: the integrator takes the site's surface "
              "and the per-patch store types; the rest is the driver's (scripts/lint/check_integrator_boundary.py)")
        sys.exit(1)
    print("check_integrator_boundary: the integrator sees one patch's records only")
