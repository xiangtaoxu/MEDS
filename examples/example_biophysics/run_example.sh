#!/usr/bin/env bash
#
# MEDS example_biophysics -- run both stages, then build the figure.
#
#   ./run_example.sh              # spin-up (if needed) + July + plot
#   ./run_example.sh --replot     # skip the model, just rebuild the figures from existing output
#
# Stage 1 is skipped automatically when its state file is already present, so re-running to tweak
# a figure costs seconds rather than the ~9 min (4-thread) / ~25 min (serial) spin-up.

set -euo pipefail
cd "$(dirname "$0")"

MEDS_BIN="${MEDS_BIN:-../../build-ifx/meds_main}"
STATE="spinup-S-20740701000000.nc"
# The plotting stack is the `meds` conda environment (environment.yml). Prefer it explicitly:
# a bare `python3` picks up whatever is first on PATH, which on a machine with several conda envs
# is routinely one with a mismatched numpy/matplotlib pair -- an ImportError at the very last step
# of a 9-minute run. Override with PYTHON=... if your environment lives elsewhere.
if [[ -z "${PYTHON:-}" ]]; then
   for cand in "$HOME/miniforge3/envs/meds/bin/python" "$CONDA_PREFIX/bin/python" python3; do
      if [[ -x "$cand" ]] && "$cand" -c "import numpy, netCDF4, matplotlib" 2>/dev/null; then
         PYTHON="$cand" ; break
      fi
   done
fi
PYTHON="${PYTHON:-python3}"

if [[ "${1:-}" != "--replot" ]]; then
   if [[ ! -f "$MEDS_BIN" ]]; then
      echo "error: MEDS binary not found at $MEDS_BIN" >&2
      echo "  build it first, e.g." >&2
      echo "    source /opt/intel/oneapi/setvars.sh" >&2
      echo "    cmake -S ../.. -B ../../build-ifx -DCMAKE_Fortran_COMPILER=ifx \\" >&2
      echo "          -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=\$HOME/miniforge3/envs/common" >&2
      echo "    cmake --build ../../build-ifx -j" >&2
      echo "  or point MEDS_BIN at an existing build." >&2
      exit 1
   fi
   if [[ ! -f "../../data/forcing/ithaca_forcing.nc" ]]; then
      echo "error: forcing file ../../data/forcing/ithaca_forcing.nc not found." >&2
      echo "  It is not tracked in git (NetCDF files are ignored). Build it with:" >&2
      echo "    python ../../scripts/download_era5land.py     # needs a CDS API key" >&2
      echo "    python ../../scripts/prep_era5land_forcing.py" >&2
      exit 1
   fi

   if [[ -f "$STATE" ]]; then
      echo "==> stage 1: $STATE already present -- skipping the 50-year spin-up"
      echo "    (delete it to force a re-run)"
   else
      echo "==> stage 1: 50-year spin-up from bare ground (2024-07-01 -> 2074-07-01)"
      echo "    no diagnostic output; the only product is $STATE.  ~9 min on 4 threads, ~25 min serial."
      "$MEDS_BIN" meds_config_spinup.toml
   fi

   echo
   echo "==> stage 2: July 2074 at hourly resolution (restart from the spin-up state)"
   mkdir -p out
   "$MEDS_BIN" meds_config_july.toml
fi

echo
echo "==> building the figures"
"$PYTHON" plot_biophysics.py     # energy:  the four temperatures
"$PYTHON" plot_carbon.py         # carbon:  GPP / NPP / Reco / NEE + canopy-air CO2
"$PYTHON" plot_soil.py           # water:   the soil moisture field, depth x time
